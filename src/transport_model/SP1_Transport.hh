#ifndef SP1_Transport_hh
#define SP1_Transport_hh

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
    using std::string;
    using std::unique_ptr;
    using std::vector;
    using namespace data_ns;
    using namespace mesh_ns;
    
    class SP1_Transport
    {
    private:
        
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

        inline double d(unsigned cell, unsigned from_group, unsigned to_group)
        {
            return d_[from_group + data_.number_of_groups() * (to_group + data_.number_of_groups() * cell)];
        }
        
    public:
        
        SP1_Transport(Data &data,
                      Mesh &mesh,
                      string problem_type = "forward");
        
        int solve();
        
        void print_scalar_flux()
        {
            std::cout << "SP1_Transport" << std::endl;
            // std::cout << "MATRIX" << std::endl;
            // std::cout << *matrix_ << std::endl;
            
            // std::cout << "RHS" << std::endl;
            // std::cout << *rhs_ << std::endl;
            
            std::cout << "LHS" << std::endl;
            std::cout << *lhs_ << std::endl;

            for (unsigned i = 0; i < number_of_edges_; ++i)
            {
                for (unsigned g = 0; g < mesh_.number_of_cells(); ++g)
                {
                    unsigned k = i + number_of_edges_ * g;

                    std::cout << (*lhs_)[k] << std::endl;
                }
            }

        }
    };
}

#endif
