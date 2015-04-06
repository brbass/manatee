#include "SPn_Transport.hh"

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <mpi.h>

#include <Amesos.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_FECrsMatrix.h>
#include <Epetra_LinearProblem.h>
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseSolver.h>
#include <Epetra_Vector.h>

#include "Data.hh"
#include "Mesh.hh"

namespace transport_ns
{
    using namespace std;
    using namespace data_ns;
    using namespace mesh_ns;
    
    SPn_Transport::
    SPn_Transport(unsigned number_of_even_moments,
                  Data &data,
                  Mesh &mesh,
                  string problem_type):
        number_of_even_moments_(number_of_even_moments),
        data_(data),
        mesh_(mesh),
        problem_type_(problem_type)
    {
        data_.check();
        mesh_.check();

        number_of_edges_ = mesh_.number_of_cells() + 1;
        
        num_global_elements_ = number_of_edges_ * number_of_even_moments_ * data_.number_of_groups();

        comm_ = unique_ptr<Epetra_MpiComm> (new Epetra_MpiComm(MPI_COMM_WORLD));
        map_ = unique_ptr<Epetra_Map> (new Epetra_Map(num_global_elements_, index_base_, *comm_));

        initialize_d();
        initialize_matrix();
        initialize_lhs();
        initialize_rhs();
        initialize_problem();
        initialize_solver();
    }

    SPn_Transport::
    ~SPn_Transport()
    {
    }
    
    int SPn_Transport::
    solve()
    {
        solver_->SymbolicFactorization();
        
        solver_->NumericFactorization();
        
        solver_->Solve();

        return 0;
    }

    int SPn_Transport::
    initialize_matrix()
    {
        num_my_elements_ = map_->NumMyElements();
        my_global_elements_.assign(map_->MyGlobalElements(), map_->MyGlobalElements() + num_my_elements_);
        
        num_entries_per_row_.resize(num_my_elements_, 0);
        
        for (unsigned i = 0; i < number_of_edges_; ++i)
        {
            num_entries_per_row_[i + 0] = 6 * data_.number_of_groups();
            
            for (unsigned m = 1; m < number_of_even_moments_ - 1; ++m)
            {
                num_entries_per_row_[i + number_of_edges_ * m] = 9 * data_.number_of_groups();
            }

            num_entries_per_row_[i + number_of_edges_ * (number_of_even_moments_ - 1)] = 6 * data_.number_of_groups();
        }
        
        matrix_ = unique_ptr<Epetra_FECrsMatrix> (new Epetra_FECrsMatrix (Copy, *map_, &num_entries_per_row_[0], true));
        
        for (unsigned i = 0; i < mesh_.number_of_cells(); ++i)
        {
            unsigned dims = mesh_.number_of_nodes() * data_.number_of_groups() * number_of_even_moments_;
            
            Epetra_SerialDenseMatrix fill_matrix (dims, dims);
            Epetra_IntSerialDenseVector fill_vector (dims);
        
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                for (unsigned m = 0; m < number_of_even_moments_; ++m)
                {
                    for (unsigned o1 = 0; o1 < mesh_.number_of_nodes(); ++o1)
                    {
                        unsigned k1 = o1 + mesh_.number_of_nodes() * (g + data_.number_of_groups() * m);
                            
                        for (unsigned o2 = 0; o2< mesh_.number_of_nodes(); ++o2)
                        {
                            unsigned k2 = o2 + mesh_.number_of_nodes() * (g + data_.number_of_groups() * m);
                                
                            fill_matrix(k1, k2) += data_.sigma_t(i, g) * mesh_.stiffness(i, o1 + 1, o2 + 1); // 1
                        }
                    }
                }
            }
            
            for (unsigned gt = 0; gt < data_.number_of_groups(); ++gt)
            {
                for (unsigned gf = 0; gf < data_.number_of_groups(); ++gf)
                {
                    for (unsigned m = 0; m < number_of_even_moments_; ++m)
                    {
                        unsigned n = 2 * m;
                        
                        double l = compute_l(i, gf, gt, n);
                        
                        for (unsigned o1 = 0; o1 < mesh_.number_of_nodes(); ++o1)
                        {
                            unsigned k1 = o1 + mesh_.number_of_nodes() * (gt + data_.number_of_groups() * m);
                            
                            for (unsigned o2 = 0; o2< mesh_.number_of_nodes(); ++o2)
                            {
                                unsigned k2 = o2 + mesh_.number_of_nodes() * (gf + data_.number_of_groups() * m);
                                
                                fill_matrix(k1, k2) += l * mesh_.stiffness(i, o1 + 1, o2 + 1, 1, 1); // 2
                            }
                            
                            fill_vector(k1) = i + o1 + number_of_edges_ * (gt + data_.number_of_groups() * m);
                        }
                        
                        if (m > 0)
                        {
                            double ll = compute_ll(i, gf, gt, n);
                            
                            for (unsigned o1 = 0; o1 < 2; ++o1)
                            {
                                unsigned k1 = o1 + mesh_.number_of_nodes() * (gt + data_.number_of_groups() * (m - 1));
                                
                                for (unsigned o2 = 0; o2< 2; ++o2)
                                {
                                    unsigned k2 = o2 + mesh_.number_of_nodes() * (gf + data_.number_of_groups() * m);
                                    fill_matrix(k1, k2) -= ll * mesh_.stiffness(i, o1 + 1, o2 + 1, 1, 1); // 4
                                }
                            }
                        }
                    
                        if (m < number_of_even_moments_ - 1)
                        {
                            double lu = compute_lu(i, gf, gt, n);
                        
                            for (unsigned o1 = 0; o1 < 2; ++o1)
                            {
                                unsigned k1 = o1 + mesh_.number_of_nodes() * (gt + data_.number_of_groups() * (m + 1));

                                for (unsigned o2 = 0; o2< 2; ++o2)
                                {
                                    unsigned k2 = o2 + mesh_.number_of_nodes() * (gf + data_.number_of_groups() * m);
                                    
                                    fill_matrix(k1, k2) -= lu * mesh_.stiffness(i, o1 + 1, o2 + 1, 1, 1); // 8
                                }
                            }
                        }
                    }
                }
            }
            
            for (unsigned gt = 0; gt < data_.number_of_groups(); ++gt)
            {
                for (unsigned gf = 0; gf < data_.number_of_groups(); ++gf)
                {
                    for (unsigned m = 0; m < number_of_even_moments_; ++m)
                    {
                        unsigned n = 2 * m;
                        
                        for (unsigned o1 = 0; o1 < mesh_.number_of_nodes(); ++o1)
                        {
                            unsigned k1 = o1 + mesh_.number_of_nodes() * (gt + data_.number_of_groups() * m);
                            
                            for (unsigned o2 = 0; o2< mesh_.number_of_nodes(); ++o2)
                            {
                                unsigned k2 = o2 + mesh_.number_of_nodes() * (gf + data_.number_of_groups() * m);
                                
                                fill_matrix(k1, k2) -= data_.sigma_s(i, gf, gt, n) * mesh_.stiffness(i, o1 + 1, o2 + 1); // 16
                            }
                        }
                    }
                }
            }
            
            matrix_->InsertGlobalValues(fill_vector, fill_matrix);
        }
        
        matrix_->GlobalAssemble();
        
        return 0;
    }
    
    int SPn_Transport::
    initialize_lhs()
    {
        lhs_ = unique_ptr<Epetra_Vector> (new Epetra_Vector(*map_));
        lhs_->PutScalar(7.7);
        
        return 0;
    }
    
    int SPn_Transport::
    initialize_rhs()
    {
        rhs_ = unique_ptr<Epetra_Vector> (new Epetra_Vector(*map_));
        rhs_->PutScalar(0.0);
        
        for (unsigned i = 0; i < mesh_.number_of_cells(); ++i)
        {
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                for (unsigned o = 0; o < mesh_.number_of_nodes(); ++o)
                {
                    unsigned k = o + i + number_of_edges_ * g;
                    
                    (*rhs_)[k] += data_.internal_source(i, g);
                }
            }
        }
        
        return 0;
    }

    int SPn_Transport::
    initialize_problem()
    {
        problem_ = unique_ptr<Epetra_LinearProblem> (new Epetra_LinearProblem(matrix_.get(), lhs_.get(), rhs_.get()));
        
        return 0;
    }

    int SPn_Transport::
    initialize_solver()
    {
        solver_ = unique_ptr<Amesos_BaseSolver> (factory_.Create(solver_type_, *problem_));
        
        if (solver_ == NULL)
        {
            cerr << "Specified solver \"" << solver_type_ << "\" is not available." << endl;
        }
        
        solver_->SetParameters(list_);
        
        return 0;
    }

    // collision and streaming
    double SPn_Transport::
    compute_l(unsigned cell, unsigned from_group, unsigned to_group, unsigned moment)
    {
        if (moment >= 2)
        {
            unsigned k = pow(moment / (moment - 1), 2) * (2 * moment - 3)/(2 * moment + 1);
            
            return pow(2 / mesh_.cell_length(cell), 2) * (k * compute_d(cell, from_group, to_group, moment - 2) + compute_d(cell, from_group, to_group, moment));
        }
        else
        {
            return pow(2 / mesh_.cell_length(cell), 2) * compute_d(cell, from_group, to_group, moment);
        }            
    }
    
    // collision and streaming, n-2
    double SPn_Transport::
    compute_ll(unsigned cell, unsigned from_group, unsigned to_group, unsigned moment)
    {
        if (moment >= 2)
        {
            unsigned k = moment / (2 * moment + 1) * (2 * moment - 3) / (moment - 1);
            
            return pow(2 / mesh_.cell_length(cell), 2) * (k * compute_d(cell, from_group, to_group, moment - 2) + compute_d(cell, from_group, to_group, moment));
        }
        else
        {
            return 0;
        }            
    }

    // collision and streaming, n+2
    double SPn_Transport::
    compute_lu(unsigned cell, unsigned from_group, unsigned to_group, unsigned moment)
    {
        unsigned k = (moment + 2) / (moment + 1);
        
        return pow(2 / mesh_.cell_length(cell), 2) * (k * compute_d(cell, from_group, to_group, moment));
    }

    int SPn_Transport::
    initialize_d()
    {
        Epetra_SerialDenseSolver solver;
        
        d_.resize(pow(data_.number_of_groups(), 2) * mesh_.number_of_cells() * (data_.number_of_scattering_moments() - 1), 0);
        
        for (unsigned i = 0; i < mesh_.number_of_cells(); ++i)
        {
            for (unsigned n = 0; n < data_.number_of_scattering_moments() - 1; ++n)
            {
                Epetra_SerialDenseMatrix invert_matrix (data_.number_of_groups(), data_.number_of_groups());
                
                for (unsigned g = 0; g < data_.number_of_groups(); ++g)
                {
                    invert_matrix(g, g) += data_.sigma_t(i, g);
                }
                
                for (unsigned gf = 0; gf < data_.number_of_groups(); ++gf)
                {
                    for (unsigned gt = 0; gt < data_.number_of_groups(); ++gt)
                    {
                        invert_matrix(gf, gt) -= data_.sigma_s(i, gf, gt, n+1);
                    }
                }

                solver.SetMatrix(invert_matrix);
                solver.Invert();
                
                for (unsigned gf = 0; gf < data_.number_of_groups(); ++gf)
                {
                    for (unsigned gt = 0; gt < data_.number_of_groups(); ++gt)
                    {
                        unsigned k = gf + data_.number_of_groups() * (gt + data_.number_of_groups() * (i + mesh_.number_of_cells() * n));

                        d_[k] = pow(n+1, 2) / ((2*n+1) * (2*n+3)) * invert_matrix(gf, gt);
                    }
                }
            }
        }

        return 0;
    }
}
