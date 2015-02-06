#include "SPn_Transport.hh"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <mpi.h>

#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
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
            
            lhs_old_[k] = lhs_[k];
        }
        
        return checksum;
    }

    int SPn_Transport::
    get_row(vector<int> &indeces, vector<double> &values, unsigned local_index, unsigned group)
    {
        unsigned i = local_to_cell_edge(local_index);
        unsigned m = local_to_moment(local_index);

        // unsigned m0 = m;
        // unsigned mN = m;
        
        // if (m > 0)
        // {
        //     m0 = m - 1;
        // }
        // if (m < number_of_even_moments_ - 1)
        // {
        //     mN = m + 1;
        // }
        
        // unsigned i0 = i;
        // unsigned iN = i;
        
        // if (i > 0)
        // {
        //     i0 = i - 1;
        // }
        // if (i < number_of_cell_edges_ - 1)
        // {
        //     i0 = i + 1;
        // }

        unsigned i1 = i;
        unsigned m1 = m;
        
        if (m != 0)
        {
            if (i != 0)
            {
                // m-1, i-1
                
                i1 = i - 1;
                m1 = m - 1;
                
                indeces.push_back(i1 + number_of_cell_edges_ * m1);
                values.push_back();
            }
            // m-1, i
            if (i != number_of_cell_edges_ - 1)
            {
                // m-1, i+1
            }
        }
        
        if (i != 0)
        {
            // m, i-1
        }
        // m, i
        if (i != number_of_cell_edges_ - 1)
        {
            // m, i+1
        }
        
        if (m != number_of_even_moments_ - 1)
        {
            if (i != 0)
            {
                // m+1, i-1
            }
            // m+1, i
            if (i != number_of_cell_edges_ - 1)
            {
                // m+1, i+1
            }
        }
    }
    
    int SPn_Transport::
    initialize_matrices()
    {
        matrix_.clear();
        row_map_.clear();
        
        num_my_elements_ = map_->NumMyElements();
        my_global_elements_.assign(map_->MyGlobalElements(), map_->MyGlobalElements() + num_my_elements_);
        
        unsigned num_diagonals = 3;
        unsigned num_off_diagonals = (num_diagonals - 1) / 2;
        
        for (int k1 = -number_of_edges_; k1 < number_of_edges_ + 1; k1 += number_of_edges_)
        {
            for (int k2 = -num_off_diagonals; k2 < num_off_diagonals; ++k2)
            {
                row_map_.push_back(k1 + k2);
            }
        }

        num_entries_per_row_.resize(num_my_elements_, 0);
        
        for (unsigned e = 0; e < num_my_elements_; ++e)
        {
            unsigned i = local_to_cell_edge(e);
            unsigned m = local_to_moment(e);

            for (int k1 = -number_of_edges_; k1 < number_of_edges_ + 1; k1 += number_of_edges_)
            {
                for (int k2 = -num_off_diagonals; k2 < num_off_diagonals; ++k2)
                {
                    if (i > 0)
                    k1 + k2;
                }
            }
            
            for (unsigned k = 0; k < row_map_.size(); ++k)
            {
                int l = i + row_map_[k];
                
                if (l >= 0 && l < number_of_edges_)
                {
                    num_entries_per_row_[i] += 1;
                }
            }
        }

        for (unsigned g = 0; g < data_.number_of_groups(); ++g)
        {
            Epetra_CrsMatrix matrix (Copy, *map_, &num_entries_per_row_[0], true);
            
            for (unsigned e = 0; i < num_my_elements_; ++i)
            {
                vector<int> indices(num_entries_per_row_[i], 0);
                vector<double> values(num_entries_per_row_[i], 0);

                if (m ==0)
                {
                    
                }
                else if (m == number_of_even_moments_ - 1)
                {
                }
                
                if (j == 0)
                {
                    values[0] = -compute_k1(j,g,n) * mesh_.stiffness_moment(j,1,1) + data_.sigma_t(j, g) * mesh_.stiffness(j,1,1);
                    values[1] = -compute_k1(j,g,n) * mesh_.stiffness_moment(j,2,1) + data_.sigma_t(j, g) * mesh_.stiffness(j,2,1);
                }
                else if (j == num_global_elements_ - 1)
                {
                    values[0] = -compute_k1(j-1,g,n) * mesh_.stiffness_moment(j-1,1,2) + data_.sigma_t(j-1, g) * mesh_.stiffness(j-1,1,2);
                    values[1] = -compute_k1(j-1,g,n) * mesh_.stiffness_moment(j-1,2,2) + data_.sigma_t(j-1, g) * mesh_.stiffness(j-1,2,2);
                }
                else
                {
                    values[0] = -compute_k1(j-1,g,n) * mesh_.stiffness_moment(j-1,1,2) + data_.sigma_t(j-1, g) * mesh_.stiffness(j-1,1,2);
                    values[1] = -compute_k1(j-1,g,n) * mesh_.stiffness_moment(j-1,2,2) + data_.sigma_t(j-1, g) * mesh_.stiffness(j-1,2,2) -compute_k1(j,g,n) * mesh_.stiffness_moment(j,1,1) + data_.sigma_t(j, g) * mesh_.stiffness(j,1,1);
                    values[2] = -compute_k1(j,g,n) * mesh_.stiffness_moment(j,2,1) + data_.sigma_t(j, g) * mesh_.stiffness(j,2,1);
                }
                matrix.InsertGlobalValues(j, num_entries_per_row_[i], &values[0], &indices[0]);
            } // e
                
            matrix.FillComplete();
            lower_matrix.FillComplete();
            upper_matrix.FillComplete();
                
            matrix_.push_back(matrix);
            lower_matrix_.push_back(lower_matrix);
            upper_matrix_.push_back(upper_matrix);
                
            for (unsigned g1 = 0; g1 < data_.number_of_groups(); ++g1)
            {
                Epetra_CrsMatrix scattering_matrix (Copy, *map_, &num_entries_per_row_[0], true);
                    
                for (unsigned i = 0; i < num_my_elements_; ++i)
                {
                    vector<int> indices(num_entries_per_row_[i], 0);
                    vector<double> values(num_entries_per_row_[i], 0);
                        
                    unsigned j = my_global_elements_[i];
                    unsigned k = 0;
                        
                    for (unsigned l = 0; l < row_map_.size(); ++l)
                    {
                        int index = j + row_map_[l];
                            
                        if (index >= 0 && index < num_global_elements_)
                        {
                            indices[k] = index;
                            k += 1;
                        }
                    }
                        
                    if (j == 0)
                    {
                        values[0] = data_.sigma_s(j, g1, g, m) * mesh_.stiffness(j,1,1);
                        values[1] = data_.sigma_s(j, g1, g, m) * mesh_.stiffness(j,2,1);
                    }
                    else if (j == num_global_elements_ - 1)
                    {
                        values[0] = data_.sigma_s(j-1, g1, g, m) * mesh_.stiffness(j-1,1,2);
                        values[1] = data_.sigma_s(j-1, g1, g, m) * mesh_.stiffness(j-1,2,2);
                    }
                    else
                    {
                        values[0] = data_.sigma_s(j-1, g1, g, m) * mesh_.stiffness(j-1,1,2);
                        values[1] = data_.sigma_s(j-1, g1, g, m) * mesh_.stiffness(j-1,2,2) + data_.sigma_s(j, g1, g, m) * mesh_.stiffness(j,1,1);
                        values[2] = data_.sigma_s(j, g1, g, m) * mesh_.stiffness(j,2,1);
                    }
                        
                    scattering_matrix.InsertGlobalValues(j, num_entries_per_row_[i], &values[0], &indices[0]);
                } // i
                
                scattering_matrix.FillComplete();
                
                scattering_matrix_.push_back(scattering_matrix);
            } // g1
        } // g
        
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
            
            unsigned j0 = my_global_elements_.front();
            unsigned i0 = 0;
            
            if (j0 == 0)
            {
                source_[g][j0] = (mesh_.stiffness(j0,1,1) + mesh_.stiffness(j0,2,1)) * data_.internal_source(j0, g);
                
                i0 = 1;
            }
            
            unsigned jn = my_global_elements_.back();
            unsigned in = num_my_elements_;
                
            if (jn == num_global_elements_ - 1)
            {
                source_[g][jn] = (mesh_.stiffness(jn-1,1,2) + mesh_.stiffness(jn-1,2,2)) * data_.internal_source(jn-1, g);
                
                in -= 1;
            }
            
            for (unsigned i = i0; i < in; ++i)
            {
                unsigned j = my_global_elements_[i];
                
                source_[g][j] = (mesh_.stiffness(j-1,1,2) + mesh_.stiffness(j-1,2,2)) * data_.internal_source(j-1, g) + (mesh_.stiffness(j,1,1) + mesh_.stiffness(j,2,1)) * data_.internal_source(j, g);
            }
        }
        
        return 0;
    }

    
    int SPn_Transport::
    calculate_rhs()
    {
        for (unsigned g = 0; g < data_.number_of_groups(); ++g)
        {
            // need to use "replace values"
            rhs_[k] = source_[g];
        }
        
        Epetra_Vector temp_vector(*map_);
        
        for (unsigned g = 0; g < data_.number_of_groups(); ++g)
        {
            for (unsigned g1 = 0; g1 < data_.number_of_groups(); ++g1)
            {
                unsigned k1 = g1 + data_.number_of_groups() *  (g + data_.number_of_groups() * m);
                    
                scattering_matrix_[g].Multiply(false, lhs_[g], temp_vector);
                    
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
    
    double SPn_Transport::
    compute_k1(unsigned cell, unsigned group, unsigned moment)
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
    
    double SPn_Transport::
    compute_k2(unsigned cell, unsigned group, unsigned moment)
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

    double SPn_Transport::
    compute_k3(unsigned cell, unsigned group, unsigned moment)
    {
        unsigned k = (moment + 2) / (moment + 1);

        return pow(2 / mesh_.cell_length(cell), 2) * (k * data_.d(cell, group, moment));
    }
}
