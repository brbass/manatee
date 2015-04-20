#ifndef SP1_Transport_hh
#define SP1_Transport_hh

#include <iomanip>
#include <iostream>
#include <cmath>
#include <memory>
#include <vector>

#include <Epetra_SerialComm.h>
//#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_FECrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_LinearProblem.h>
#include <Amesos.h>

#include "Data.hh"
#include "Mesh.hh"
#include "Transport_Model.hh"

namespace transport_ns
{
    using std::string;
    using std::unique_ptr;
    using std::vector;
    using namespace data_ns;
    using namespace mesh_ns;
    
    class SP1_Transport: public Transport_Model
    {
    private:
        
        unsigned number_of_edges_;
        unsigned max_num_iterations_ = 1000;
        
        double tolerance = 1e-8;
        
        Data data_;
        Mesh mesh_;
        
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
        unique_ptr<Epetra_SerialComm> comm_;
        //unique_ptr<Epetra_MpiComm> comm_;
        unique_ptr<Epetra_FECrsMatrix> matrix_;
        unique_ptr<Epetra_Vector> lhs_;
        unique_ptr<Epetra_Vector> rhs_;
        unique_ptr<Epetra_LinearProblem> problem_;
        
        vector<double> d_;
        vector<double> phi_;
        
        string solver_type_ = "Klu";
        unique_ptr<Amesos_BaseSolver> solver_;
        Amesos factory_;
        
        Teuchos::ParameterList list_;
        
        void initialize_d();
        void initialize_matrix();
        void initialize_lhs();
        void initialize_rhs();
        void initialize_problem();
        void initialize_solver();
        void initialize_transport();

        void transfer_scalar_flux();

        inline double d(unsigned cell, unsigned from_group, unsigned to_group)
        {
            return d_[from_group + data_.number_of_groups() * (to_group + data_.number_of_groups() * cell)];
        }
        
    public:
        
        SP1_Transport(unsigned &number_of_cells,
                      unsigned &number_of_groups,
                      unsigned &number_of_scattering_moments,
                      double &side_length,
                      vector<double> &internal_source,
                      vector<double> &boundary_sources,
                      vector<double> &sigma_t,
                      vector<double> &sigma_s,
                      vector<double> &nu_sigma_f,
                      vector<double> &chi,
                      vector<string> &boundary_conditions,
                      string problem_type = "forward");
        
        virtual void solve();
        
        virtual void print_scalar_flux()
        {
            using namespace std;

            const int w = 8;
            
            cout << "SP1_Transport" << endl;
            // cout << "MATRIX" << endl;
            // cout << *matrix_ << endl;
            
            // cout << "RHS" << endl;
            // cout << *rhs_ << endl;
            
            // cout << "LHS" << endl;
            // cout << *lhs_ << endl;

            cout << left;
            cout << setw(w) << "cell" << setw(w) << "node" << setw(w) << "group" << setw(w) << "phi" <<  endl;
            
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                for (unsigned i = 0; i < mesh_.number_of_cells(); ++i)
                {
                    for (unsigned n = 0; n < mesh_.number_of_nodes(); ++n)
                    {
                        unsigned k = i + n + number_of_edges_ * g;
                        
                        cout << setw(w) << i << setw(w) << n << setw(w) << g << setw(w) << (*lhs_)[k] << endl;
                    }
                }
            }

            cout << endl;
        }

        vector<double> get_scalar_flux()
        {
            return phi_;
        }
    };
}

#endif
