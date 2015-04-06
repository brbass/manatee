#ifndef SPn_Transport_hh
#define SPn_Transport_hh

#include <iostream>
#include <cmath>
#include <memory>
#include <vector>

#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_FECrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_LinearProblem.h>
#include <Amesos.h>

#include "Data.hh"
#include "Mesh.hh"

namespace transport_ns
{
    using std::unique_ptr;
    using std::vector;
    using namespace data_ns;
    using namespace mesh_ns;
    
    class SPn_Transport
    {
    private:
        
        unsigned number_of_even_moments_;
        unsigned number_of_edges_;
        unsigned max_num_iterations_ = 1000;
        
        double tolerance = 1e-8;
        
        Data &data_;
        Mesh &mesh_;

        string problem_type_;
        
        vector<double> solution_;
        
        inline double psi(unsigned cell, unsigned group, unsigned moment);
        
        const int index_base_ = 0;
        int num_global_elements_ = 0; // number of elements (rows) on all processors
        int num_my_elements_ = 0; // number of elements on processor
        
        vector<int> num_entries_per_row_;
        vector<int> my_global_elements_; // global indeces of elements (rows)
        vector<int> row_map_; // position of elements from diagonal
        
        unique_ptr<Epetra_Map> map_;
        unique_ptr<Epetra_MpiComm> comm_;
        unique_ptr<Epetra_FECrsMatrix> matrix_;
        unique_ptr<Epetra_Vector> lhs_;
        unique_ptr<Epetra_Vector> rhs_;
        unique_ptr<Epetra_LinearProblem> problem_;
        
        vector<double> d_;
        
        string solver_type_ = "Klu";
        unique_ptr<Amesos_BaseSolver> solver_;
        Amesos factory_;

        Teuchos::ParameterList list_;
        
        int initialize_d();
        int initialize_matrix();
        int initialize_lhs();
        int initialize_rhs();
        int initialize_problem();
        int initialize_solver();
        int initialize_transport();
        
        inline int local_to_cell_edge(int local_index)
        {
            return my_global_elements_[local_index] % number_of_edges_;
        }

        inline int local_to_moment(int local_index)
        {
            return my_global_elements_[local_index] / number_of_edges_;
        }

        inline double compute_d(unsigned cell, unsigned from_group, unsigned to_group, unsigned moment)
        {
            return d_[from_group + data_.number_of_groups() * (to_group + data_.number_of_groups() * (cell + mesh_.number_of_cells() * moment))];
        }
        
        double compute_l(unsigned cell, unsigned from_group, unsigned to_group, unsigned moment);
        double compute_ll(unsigned cell, unsigned from_group, unsigned to_group,unsigned moment);
        double compute_lu(unsigned cell, unsigned from_group, unsigned to_group, unsigned moment);
        
    public:

        SPn_Transport(unsigned number_of_even_moments,
                      Data &data,
                      Mesh &mesh,
                      string problem_type = "forward");

        ~SPn_Transport();
        
        int solve();

        void print_scalar_flux()
        {
            std::cout << "SPn_Transport" << std::endl;
            std::cout << "MATRIX" << std::endl;
            std::cout << *matrix_ << std::endl;
            
            std::cout << "RHS" << std::endl;
            std::cout << *rhs_ << std::endl;
            
            std::cout << "LHS" << std::endl;
            std::cout << *lhs_ << std::endl;
            
            // plot_scalar_flux();
        }

        void plot_scalar_flux()
        {
            unsigned lines = 20;
            
            if (number_of_edges_ < 500)
            {
                for (unsigned g = 0; g < data_.number_of_groups(); ++g)
                {
                    std::cout << "Group " << g << " Normalized Scalar Flux ";
                    
                    vector<double> temp1(number_of_edges_, 0);
                    
                    double max = 0;
                    for (unsigned i = 0; i < number_of_edges_; ++i)
                    {
                        unsigned k = 0 + i + number_of_edges_ * (g + 0);
                        
                        temp1[i] = (*lhs_)[k];
                        
                        if (temp1[i] > max)
                        {
                            max = temp1[i];
                        }
                    }

                    std::cout << "(phi_max = " << max << ")" << std::endl;
                    
                    vector<unsigned> temp2(number_of_edges_,0);
                    
                    for (unsigned i = 0; i < number_of_edges_; ++i)
                    {
                        temp1[i] /= max;
                        
                        temp2[i] = std::ceil(temp1[i] * lines);
                    }
                    
                    
                    for (int l = lines; l >= 0; --l)
                    {
                        for (unsigned i = 0; i < number_of_edges_; ++i)
                        {
                            if (temp2[i] == l)
                            {
                                std::cout << ".";
                            }
                            else
                            {
                                std::cout << " ";
                            }
                        }
                        std::cout << std::endl;
                    }

                    std::cout << std::endl;
                }
            }
        }
    };
}
#endif
