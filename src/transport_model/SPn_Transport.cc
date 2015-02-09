#include "SPn_Transport.hh"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <mpi.h>

#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_FECrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_LinearProblem.h>
#include <Amesos.h>

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
                  Mesh &mesh):
        number_of_even_moments_(number_of_even_moments),
        data_(data),
        mesh_(mesh)
    {
        data_.check();
        mesh_.check();

        number_of_edges_ = mesh_.number_of_cells() + 1;
        
        num_global_elements_ = number_of_edges_ * number_of_even_moments_;
        
        comm_ = new Epetra_MpiComm(MPI_COMM_WORLD);
        map_ = new Epetra_Map(num_global_elements_, index_base_, *comm_);

        initialize_matrices();
        initialize_lhs();
        initialize_rhs();
        initialize_source();
        initialize_problem();
        initialize_solver();
    }

    SPn_Transport::
    ~SPn_Transport()
    {
        delete map_;
        delete comm_;
        for (unsigned i = 0; i < solver_.size(); ++i)
        {
            delete solver_[i];
        }
    }
    
    int SPn_Transport::
    solve()
    {
        unsigned iterations = 0;
        
        for (unsigned it = 0; it < max_num_iterations_; ++ it)
        {
            calculate_rhs();
            
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                solver_[g]->Solve();
            }
            
            if (!check_convergence())
            {
                iterations = it + 1;
                
                break;
            }
        }

        if (iterations == 0)
        {
            iterations = max_num_iterations_;
        }
        
        cout << "iterations: " << iterations << endl;
        
        return 0;
    }

    int SPn_Transport::
    check_convergence()
    {
        int checksum = 0;
        
        Epetra_Vector temp_vector(*map_);

        for (unsigned g = 0; g < data_.number_of_groups(); ++g)
        {
            double lhs_norm = 0;
            
            temp_vector = lhs_[g];
            temp_vector.Scale(-1);
            temp_vector.SumIntoGlobalValues(num_my_elements_, &lhs_old_[g][my_global_elements_[0]], &my_global_elements_[0]);
            temp_vector.Norm1(&lhs_norm);
            
            if (lhs_norm > tolerance)
            {
                checksum += 1;
            }
            
            lhs_old_[g] = lhs_[g];
        }
        
        return checksum;
    }

    int SPn_Transport::
    initialize_matrices()
    {
        matrix_.clear();
        row_map_.clear();
        
        num_my_elements_ = map_->NumMyElements();
        my_global_elements_.assign(map_->MyGlobalElements(), map_->MyGlobalElements() + num_my_elements_);

        num_entries_per_row_.resize(num_my_elements_, 0);
        
        for (unsigned i = 0; i < number_of_edges_; ++i)
        {
            num_entries_per_row_[i + 0] = 6;
            
            for (unsigned m = 1; m < number_of_even_moments_ - 1; ++m)
            {
                num_entries_per_row_[i + number_of_edges_ * m] = 9;
            }

            num_entries_per_row_[i + number_of_edges_ * (number_of_even_moments_ - 1)] = 6;
        }

        for (unsigned g = 0; g < data_.number_of_groups(); ++g)
        {
            Epetra_FECrsMatrix matrix (Copy, *map_, &num_entries_per_row_[0], true);
            Epetra_SerialDenseMatrix fill_matrix (2, 2);
            Epetra_IntSerialDenseVector fill_vector (2);
                    
            for (unsigned m = 0; m < number_of_even_moments_; ++m)
            {
                unsigned n = 2 * m;
                
                for (unsigned i = 0; i < mesh_.number_of_cells(); ++i)
                {
                    double l = compute_l(i, g, n);
                    
                    for (unsigned o1 = 0; o1 < 2; ++o1)
                    {
                        for (unsigned o2 = 0; o2< 2; ++o2)
                        {
                            fill_matrix(o1, o2) = data_.sigma_t(i, g) * mesh_.stiffness(i, o1 + 1, o2 + 1) - l * mesh_.stiffness_moment(i, o1 + 1, o2 + 1);
                        }
                        
                        fill_vector(o1) = i + o1 + number_of_edges_ * m;
                    }
                    
                    matrix.InsertGlobalValues(fill_vector, fill_matrix);
                    
                    if (m > 0)
                    {
                        double ll = compute_ll(i, g, n);
                        
                        for (unsigned o1 = 0; o1 < 2; ++o1)
                        {
                            for (unsigned o2 = 0; o2< 2; ++o2)
                            {
                                fill_matrix(o1, o2) = -ll * mesh_.stiffness_moment(i, o1 + 1, o2 + 1);
                            }
                            
                            fill_vector(o1) = i + o1 + number_of_edges_ * (m - 1);
                        }

                        matrix.InsertGlobalValues(fill_vector, fill_matrix);
                    }
                    
                    if (m < number_of_even_moments_ - 1)
                    {
                        double lu = compute_lu(i, g, n);
                        
                        for (unsigned o1 = 0; o1 < 2; ++o1)
                        {
                            for (unsigned o2 = 0; o2< 2; ++o2)
                            {
                                fill_matrix(o1, o2) = -lu * mesh_.stiffness_moment(i, o1 + 1, o2 + 1);
                            }

                            fill_vector(o1) = i + o1 + number_of_edges_ * (m + 1);
                        }

                        matrix.InsertGlobalValues(fill_vector, fill_matrix);
                    }
                }
            }
            
            matrix.GlobalAssemble();
            
            matrix_.push_back(matrix);
        }
        
        for (unsigned g = 0; g < data_.number_of_groups(); ++g)
        {
            for (unsigned g1 = 0; g1 < data_.number_of_groups(); ++g1)
            {
                Epetra_FECrsMatrix matrix (Copy, *map_, &num_entries_per_row_[0], true);
                Epetra_SerialDenseMatrix fill_matrix (2, 2);
                Epetra_IntSerialDenseVector fill_vector (2);
                
                for (unsigned m = 0; m < number_of_even_moments_; ++m)
                {
                    unsigned n = 2 * m;
                
                    for (unsigned i = 0; i < mesh_.number_of_cells(); ++i)
                    {
                        for (unsigned o1 = 0; o1 < 2; ++o1)
                        {
                            for (unsigned o2 = 0; o2< 2; ++o2)
                            {
                                fill_matrix(o1, o2) = data_.sigma_s(i, g, g1, n) * mesh_.stiffness(i, o1 + 1, o2 + 1);
                            }
                        
                            fill_vector(o1) = i + o1 + number_of_edges_ * m;
                        }
                    
                        matrix.InsertGlobalValues(fill_vector, fill_matrix);
                    }
                }
                
                matrix.GlobalAssemble();
                
                scattering_matrix_.push_back(matrix);
            }
        }
        
        return 0;
    }
    
    int SPn_Transport::
    initialize_lhs()
    {
        lhs_.clear();
        lhs_old_.clear();
        
        for (unsigned g = 0; g < data_.number_of_groups(); ++g)
        {
            lhs_.emplace_back(*map_);
            lhs_old_.emplace_back(*map_);
                
            lhs_.back().PutScalar(1.0);
            lhs_old_.back().PutScalar(0.0);
        }
        
        return 0;
    }
    
    int SPn_Transport::
    initialize_rhs()
    {
        rhs_.clear();
        
        for (unsigned g = 0; g < data_.number_of_groups(); ++g)
        {
            rhs_.emplace_back(*map_);
                
            rhs_.back().PutScalar(1.0);
        }

        return 0;
    }

    int SPn_Transport::
    initialize_source()
    {
        source_.clear();
        
        for (unsigned g = 0; g < data_.number_of_groups(); ++g)
        {
            source_.emplace_back(*map_);
            source_[g].PutScalar(0.0);

            for (unsigned i = 0; i < mesh_.number_of_cells(); ++i)
            {
                source_[g][i/* +0 */] += data_.internal_source(i, g);
                source_[g][i+1/* +0 */] += data_.internal_source(i, g);
            }

        }

        return 0;
    }

    
    int SPn_Transport::
    calculate_rhs()
    {
        for (unsigned g = 0; g < data_.number_of_groups(); ++g)
        {
            rhs_[g] = source_[g];
        }
        
        Epetra_Vector temp_vector(*map_);
        
        for (unsigned g = 0; g < data_.number_of_groups(); ++g)
        {
            for (unsigned g1 = 0; g1 < data_.number_of_groups(); ++g1)
            {
                unsigned k = g1 + data_.number_of_groups() *  g;
                    
                scattering_matrix_[k].Multiply(false, lhs_[g], temp_vector);
                
                rhs_[g].SumIntoGlobalValues(num_my_elements_, &temp_vector[my_global_elements_[0]], &my_global_elements_[0]);
            }
                
            problem_[g].SetRHS(&rhs_[g]);
        }

        return 0;
    }

    int SPn_Transport::
    initialize_problem()
    {
        problem_.clear();
        
        for (unsigned g = 0; g < data_.number_of_groups(); ++g)
        {
            problem_.emplace_back(&matrix_[g], &lhs_[g], &rhs_[g]);
        }

        return 0;
    }

    int SPn_Transport::
    initialize_solver()
    {
        solver_.clear();
        
        for (unsigned g = 0; g < data_.number_of_groups(); ++g)
        {
            solver_.push_back(factory_.Create(solver_type_, problem_[g]));
                
            if (solver_.back() == NULL)
            {
                cerr << "Specified solver \"" << solver_type_ << "\" is not available." << endl;
            }
        
            solver_.back()->SetParameters(list_);
                
            solver_.back()->SymbolicFactorization();
                
            solver_.back()->NumericFactorization();
        }

        return 0;
    }

    // collision and streaming
    double SPn_Transport::
    compute_l(unsigned cell, unsigned group, unsigned moment)
    {
        if (moment >= 2)
        {
            unsigned k = pow(moment / (moment - 1), 2) * (2 * moment - 3)/(2 * moment + 1);
            
            return pow(2 / mesh_.cell_length(cell), 2) * (k * data_.d(cell, group, moment - 2) + data_.d(cell, group, moment));
        }
        else
        {
            return pow(2 / mesh_.cell_length(cell), 2) * data_.d(cell, group, moment);
        }            
    }

    // collision and streaming, n-2
    double SPn_Transport::
    compute_ll(unsigned cell, unsigned group, unsigned moment)
    {
        if (moment >= 2)
        {
            unsigned k = moment / (2 * moment + 1) * (2 * moment - 3) / (moment - 1);
            
            return pow(2 / mesh_.cell_length(cell), 2) * (k * data_.d(cell, group, moment - 2) + data_.d(cell, group, moment));
        }
        else
        {
            return 0;
        }            
    }

    // collision and streaming, n+2
    double SPn_Transport::
    compute_lu(unsigned cell, unsigned group, unsigned moment)
    {
        unsigned k = (moment + 2) / (moment + 1);

        return pow(2 / mesh_.cell_length(cell), 2) * (k * data_.d(cell, group, moment));
    }
}
