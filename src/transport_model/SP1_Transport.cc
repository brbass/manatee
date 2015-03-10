#include "SP1_Transport.hh"

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
    
    SP1_Transport::
    SP1_Transport(Data &data,
                  Mesh &mesh):
        data_(data),
        mesh_(mesh)
    {
        data_.check();
        mesh_.check();

        number_of_edges_ = mesh_.number_of_cells() + 1;
        
        num_global_elements_ = number_of_edges_ * data_.number_of_groups();

        comm_ = unique_ptr<Epetra_MpiComm> (new Epetra_MpiComm(MPI_COMM_WORLD));
        map_ = unique_ptr<Epetra_Map> (new Epetra_Map(num_global_elements_, index_base_, *comm_));

        initialize_d();
        initialize_matrix();
        initialize_lhs();
        initialize_rhs();
        initialize_problem();
        initialize_solver();
    }

    SP1_Transport::
    ~SP1_Transport()
    {
    }
    
    int SP1_Transport::
    solve()
    {
        solver_->SymbolicFactorization();
        
        solver_->NumericFactorization();
        
        solver_->Solve();
        
        solver_->PrintStatus();
        
        return 0;
    }
    
    int SP1_Transport::
    initialize_matrix()
    {
        num_my_elements_ = map_->NumMyElements();
        my_global_elements_.assign(map_->MyGlobalElements(), map_->MyGlobalElements() + num_my_elements_);
        
        num_entries_per_row_.resize(num_my_elements_, 0);
        
        for (unsigned i = 0; i < number_of_edges_; ++i)
        {
            num_entries_per_row_[i + 0] = 3 * data_.number_of_groups();
        }
        
        matrix_ = unique_ptr<Epetra_FECrsMatrix> (new Epetra_FECrsMatrix (Copy, *map_, &num_entries_per_row_[0], true));
        
        unsigned n = 0;
        
        for (unsigned i = 0; i < mesh_.number_of_cells(); ++i)
        {
            unsigned dims = mesh_.number_of_nodes() * data_.number_of_groups();
            
            Epetra_SerialDenseMatrix fill_matrix (dims, dims);
            Epetra_IntSerialDenseVector fill_vector (dims);
            
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                for (unsigned o1 = 0; o1 < mesh_.number_of_nodes(); ++o1)
                {
                    unsigned k1 = o1 + mesh_.number_of_nodes() * g;
                    
                    for (unsigned o2 = 0; o2< mesh_.number_of_nodes(); ++o2)
                    {
                        unsigned k2 = o2 + mesh_.number_of_nodes() * g;
                                
                        fill_matrix(k1, k2) += data_.sigma_t(i, g) * mesh_.stiffness(i, o1 + 1, o2 + 1); // 1
                    }
                }
            }
            
            for (unsigned gt = 0; gt < data_.number_of_groups(); ++gt)
            {
                for (unsigned gf = 0; gf < data_.number_of_groups(); ++gf)
                {
                    double l = compute_l(i, gf, gt);
                    
                    for (unsigned o1 = 0; o1 < mesh_.number_of_nodes(); ++o1)
                    {
                        unsigned k1 = o1 + mesh_.number_of_nodes() * gt;
                            
                        for (unsigned o2 = 0; o2< mesh_.number_of_nodes(); ++o2)
                        {
                            unsigned k2 = o2 + mesh_.number_of_nodes() * gf;
                                
                            fill_matrix(k1, k2) -= l * mesh_.stiffness_moment(i, o1 + 1, o2 + 1); // 2
                        }
                            
                        fill_vector(k1) = i + o1 + number_of_edges_ * gt;
                    }
                }
            }

            for (unsigned gt = 0; gt < data_.number_of_groups(); ++gt)
            {
                for (unsigned gf = 0; gf < data_.number_of_groups(); ++gf)
                {
                    for (unsigned o1 = 0; o1 < mesh_.number_of_nodes(); ++o1)
                    {
                        unsigned k1 = o1 + mesh_.number_of_nodes() * gt;
                            
                        for (unsigned o2 = 0; o2< mesh_.number_of_nodes(); ++o2)
                        {
                            unsigned k2 = o2 + mesh_.number_of_nodes() * gf;
                                
                            fill_matrix(k1, k2) -= data_.sigma_s(i, gf, gt, n) * mesh_.stiffness(i, o1 + 1, o2 + 1); // 16
                        }
                    }
                }
            }
            
            // if (i == 0)
            // {
            //     for (unsigned gt = 0; gt < data_.number_of_groups(); ++gt)
            //     {
            //         unsigned o = 0;
                    
            //         unsigned k1 = o + mesh_.number_of_nodes() * gt;

            //         for (unsigned gf = 0; gf < data_.number_of_groups(); ++gf)
            //         {
            //             for (unsigned o2 = 0; o2< mesh_.number_of_nodes(); ++o2)
            //             {
            //                 unsigned k2 = o2 + mesh_.number_of_nodes() * gf;
                            
            //                 fill_matrix(k1, k2) = 0;
            //             }
            //         }

            //         fill_matrix(k1, k1) = 1;
            //     }
            // }
            // else if (i == mesh_.number_of_cells() - 1)
            // {
            //     for (unsigned gt = 0; gt < data_.number_of_groups(); ++gt)
            //     {
            //         unsigned o = 1;
                        
            //         unsigned k1 = o + mesh_.number_of_nodes() * gt;
                            
            //         for (unsigned gf = 0; gf < data_.number_of_groups(); ++gf)
            //         {
            //             for (unsigned o2 = 0; o2< mesh_.number_of_nodes(); ++o2)
            //             {
            //                 unsigned k2 = o2 + mesh_.number_of_nodes() * gf;
                            
            //                 fill_matrix(k1, k2) = 0;
            //             }
            //         }

            //         fill_matrix(k1, k1) = 1;
            //     }
            // }
            
            matrix_->InsertGlobalValues(fill_vector, fill_matrix);
        }
        
        matrix_->GlobalAssemble();
        
        return 0;
    }
    
    int SP1_Transport::
    initialize_lhs()
    {
        lhs_ = unique_ptr<Epetra_Vector> (new Epetra_Vector(*map_));
        lhs_->PutScalar(7.7);
        
        return 0;
    }
    
    int SP1_Transport::
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

                    // if ((i != 0 || o != 0) && (i != mesh_.number_of_cells() - 1 || o != 1))
                    // {
                        (*rhs_)[k] += data_.internal_source(i, g);
                    // }
                }
            }
        }

        return 0;
    }

    int SP1_Transport::
    initialize_problem()
    {
        problem_ = unique_ptr<Epetra_LinearProblem> (new Epetra_LinearProblem(matrix_.get(), lhs_.get(), rhs_.get()));
        
        return 0;
    }

    int SP1_Transport::
    initialize_solver()
    {
        solver_ = unique_ptr<Amesos_BaseSolver> (factory_.Create(solver_type_, *problem_));
        
        if (solver_ == NULL)
        {
            cerr << "Specified solver \"" << solver_type_ << "\" is not available." << endl;
        }
        
        list_.set("PrintTiming", true);
        list_.set("PrintStatus", true);
        
        solver_->SetParameters(list_);
        
        return 0;
    }

    // collision and streaming
    double SP1_Transport::
    compute_l(unsigned cell, unsigned from_group, unsigned to_group)
    {
        return pow(2 / mesh_.cell_length(cell), 2) * compute_d(cell, from_group, to_group);
    }
    
    int SP1_Transport::
    initialize_d()
    {
        unsigned n = 0;
        
        Epetra_SerialDenseSolver solver;
        
        d_.resize(pow(data_.number_of_groups(), 2) * mesh_.number_of_cells(), 0);
        
        for (unsigned i = 0; i < mesh_.number_of_cells(); ++i)
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
                    unsigned k = gf + data_.number_of_groups() * (gt + data_.number_of_groups() * i);

                    d_[k] = pow(n+1, 2) / ((2*n+1) * (2*n+3)) * invert_matrix(gf, gt);
                }
            }
        }
        
        return 0;
    }
}
