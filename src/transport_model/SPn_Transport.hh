#ifndef SPn_Transport_hh
#define SPn_Transport_hh

#include <iostream>
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
        
        vector<double> solution_;
        
        inline double psi(unsigned cell, unsigned group, unsigned moment);
        
        const int index_base_ = 0;
        int num_global_elements_ = 0; // number of elements (rows) on all processors
        int num_my_elements_ = 0; // number of elements on processor
        
        vector<int> num_entries_per_row_;
        vector<int> my_global_elements_; // global indeces of elements (rows)
        vector<int> row_map_; // position of elements from diagonal
        
        Epetra_Map *map_ = 0;
        Epetra_MpiComm *comm_ = 0;
        vector<Epetra_FECrsMatrix> matrix_;
        vector<Epetra_CrsMatrix> scattering_matrix_;
        vector<Epetra_Vector> lhs_;
        vector<Epetra_Vector> lhs_old_;
        vector<Epetra_Vector> rhs_;
        vector<Epetra_Vector> source_;
        vector<Epetra_LinearProblem> problem_;
        
        string solver_type_ = "Klu";
        
        vector<Amesos_BaseSolver*> solver_;
        Amesos factory_;
        Teuchos::ParameterList list_;
        
        int initialize_matrices();
        int initialize_lhs();
        int initialize_rhs();
        int initialize_source();
        int initialize_problem();
        int initialize_solver();
        int initialize_transport();
        int calculate_rhs();
        int check_convergence();
        
        inline int local_to_cell_edge(int local_index)
        {
            return my_global_elements_[local_index] % number_of_edges_;
        }
        inline int local_to_moment(int local_index)
        {
            return my_global_elements_[local_index] / number_of_edges_;
        }
        
        double compute_l(unsigned cell, unsigned group, unsigned moment);
        double compute_ll(unsigned cell, unsigned group, unsigned moment);
        double compute_lu(unsigned cell, unsigned group, unsigned moment);

    public:

        SPn_Transport(unsigned number_of_even_moments,
                      Data &data,
                      Mesh &mesh);

        ~SPn_Transport();
        
        int solve();

        void print_scalar_flux()
        {
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                std::cout << "Group " << g << " scalar flux" << std::endl;
                std::cout << lhs_[g] << std::endl;

                // std::cout << "Group " << g << " lhs" << std::endl;
                // std::cout << rhs_[g] << std::endl;

                // std::cout << "Group " << g << " matrix" << std::endl;
                // std::cout << rhs_[g] << std::endl;

                // std::cout << "Group " << g << " source" << std::endl;
                // std::cout << source_[g] << std::endl;

                // for (unsigned g1 = 0; g1 < data_.number_of_groups(); ++g1)
                // {
                //     std::cout << "From " << g << " to " << g1 << " group scattering" << std::endl;
                //     std::cout << scattering_matrix_[g + data_.number_of_groups() * g1] << std::endl;
                // }
            }
        }
        
    };
}
#endif
