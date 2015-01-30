#include "SPn_Transport.hh"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <mpi.h>

#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>
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
    SPn_Transport(unsigned number_of_moments,
                  Data &data,
                  Mesh &mesh):
        number_of_moments_(number_of_moments),
        data_(data),
        mesh_(mesh)
    {
        data_.check();
        mesh_.check();
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
        return 0;
    }

    int SPn_Transport::
    initialize_matrix()
    {
        matrix_.clear();
        
        num_my_elements_ = map_->NumMyElements();
        my_global_elements_.assign(map_->MyGlobalElements(), map_->MyGlobalElements() + num_my_elements_);
        
        unsigned num_diagonals = 3;
        
        row_map_.resize(num_diagonals, 0);
        row_map_[0] = -1;
        row_map_[1] = 0;
        row_map_[2] = 1;
        
        num_entries_per_row_.resize(num_my_elements_, 0);

        for (int i = 0; i < num_my_elements_; ++i)
        {
            int j = my_global_elements_[i];

            for (int k = 0; k < row_map_.size(); ++k)
            {
                if (j + row_map_[k] >= 0 && j + row_map_[k] < num_global_elements_)
                {
                    num_entries_per_row_[i] += 1;
                }
            }
        }
        
        for (unsigned m = 0; m < number_of_moments_; ++m)
        {
            unsigned n = 2 * m;
            
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                Epetra_CrsMatrix matrix (Copy, *map_, &num_entries_per_row_[0], true);
                
                for (unsigned i = 0; i < num_my_elements_; ++i)
                {
                    vector<int> indices(num_entries_per_row_[i], 0);
                    vector<double> values(num_entries_per_row_[i], 0);
                    
                    unsigned j = my_global_elements_[i];
                    unsigned k = 0;
                    
                    for (int l = 0; l < row_map_.size(); ++l)
                    {
                        int index = j + row_map_[l];
                        
                        if (index >= 0 && index < num_global_elements_)
                        {
                            indices[k] == index;
                            k += 1;
                        }
                    }
                    
                    if (j == 0)
                    {
                        double k1 = compute_k1(j, g, n);
                        
                        values[0] = k1 * mesh_.stiffness_moment(j,1,1) + data_.sigma_t(j, g) * mesh_.stiffness(j,1,1);
                        values[1] = k1 * mesh_.stiffness_moment(j,2,1) + data_.sigma_t(j, g) * mesh_.stiffness(j,2,1);
                    }
                    else if (j == num_global_elements_ - 1)
                    {
                        double k1 = compute_k1(j-1, g, n);
                        
                        values[0] = k1 * mesh_.stiffness_moment(j-1,1,2) + data_.sigma_t(j-1, g) * mesh_.stiffness(j-1,1,2);
                        values[1] = k1 * mesh_.stiffness_moment(j-1,2,2) + data_.sigma_t(j-1, g) * mesh_.stiffness(j-1,2,2);
                    }
                    else
                    {
                        double k1a = compute_k1(j-1, g, n);
                        double k1b = compute_k1(j, g, n);
                        
                        values[0] = k1a * mesh_.stiffness_moment(j-1,1,2) + data_.sigma_t(j-1, g) * mesh_.stiffness(j-1,1,2);
                        values[1] = k1a * mesh_.stiffness_moment(j-1,2,2) + data_.sigma_t(j-1, g) * mesh_.stiffness(j-1,2,2) + k1b * mesh_.stiffness_moment(j,1,1) + data_.sigma_t(j, g) * mesh_.stiffness(j,1,1);
                        values[2] = k1b * mesh_.stiffness_moment(j,2,1) + data_.sigma_t(j, g) * mesh_.stiffness(j,2,1);
                    }
                    
                    matrix.InsertGlobalValues(j, num_entries_per_row_[i], &values[0], &indices[0]);
                }

                matrix.FillComplete();

                matrix_.push_back(matrix);
            }
        }

        return 0;
    }

    int SPn_Transport::
    initialize_lhs()
    {
        lhs_.clear();
        
        for (unsigned m = 0; m < number_of_moments_; ++m)
        {
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                lhs_.emplace_back(*map_, num_vectors_);
                
                lhs_.back().PutScalar(1.0);
            }
        }

        return 0;
    }

    int SPn_Transport::
    initialize_rhs()
    {
        rhs_.clear();
        
        for (unsigned m = 0; m < number_of_moments_; ++m)
        {
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                rhs_.emplace_back(*map_, num_vectors_);
                
                rhs_.back().PutScalar(1.0);
            }
        }

        return 0;
    }

    int SPn_Transport::
    calculate_rhs()
    {
        return 0;
    }

    int SPn_Transport::
    initialize_problem()
    {
        problem_.clear();
        
        for (unsigned m = 0; m < number_of_moments_; ++m)
        {
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                unsigned k = g + data_.number_of_groups() * m;

                problem_.emplace_back(&matrix_[k], &lhs_[k], &rhs_[k]);
            }
        }

        return 0;
    }

    int SPn_Transport::
    initialize_solver()
    {
        solver_.clear();
        
        for (unsigned m = 0; m < number_of_moments_; ++m)
        {
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                unsigned k = g + data_.number_of_groups() * m;

                solver_.push_back(factory_.Create(solver_type_, problem_[k]));
                
                if (solver_.back() == NULL)
                {
                    cerr << "Specified solver \"" << solver_type_ << "\" is not available." << endl;
                }
        
                solver_.back()->SetParameters(list_);
                
                solver_.back()->SymbolicFactorization();
                
                solver_.back()->NumericFactorization();
            }
        }

        return 0;
    }
    
    int SPn_Transport::
    initialize_transport()
    {
        num_global_elements_ = mesh_.number_of_cells() + 1;
        
        comm_ = new Epetra_MpiComm(MPI_COMM_WORLD);
        map_ = new Epetra_Map(num_global_elements_, index_base_, *comm_);
        
        initialize_matrix();
        initialize_lhs();
        initialize_rhs();
        initialize_problem();
        initialize_solver();

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
