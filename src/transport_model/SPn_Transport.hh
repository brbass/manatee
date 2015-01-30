#ifndef SPn_Transport_hh
#define SPn_Transport_hh

#include <vector>

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
    using std::vector;
    using namespace data_ns;
    using namespace mesh_ns;
    
    class SPn_Transport
    {
    public:

        SPn_Transport(unsigned number_of_moments,
                      Data &data,
                      Mesh &mesh);

        ~SPn_Transport();
        
        int solve();
        
    private:
        
        unsigned number_of_moments_;
        
        Data &data_;
        Mesh &mesh_;
        
        vector<double> solution_;
        
        inline double psi(unsigned cell, unsigned group, unsigned moment);
        
        const int index_base_ = 0;
        int num_global_elements_ = 0; // number of elements (rows) on all processors
        int num_my_elements_ = 0; // number of elements on processor
        const int num_vectors_ = 1; // number of solution vectors
        
        vector<int> num_entries_per_row_;
        vector<int> my_global_elements_; // global indeces of elements (rows)
        vector<int> row_map_; // position of elements from diagonal
        
        Epetra_Map *map_ = 0;
        Epetra_MpiComm *comm_ = 0;
        vector<Epetra_CrsMatrix> matrix_;
        vector<Epetra_MultiVector> lhs_;
        vector<Epetra_MultiVector> rhs_;
        vector<Epetra_LinearProblem> problem_;
        
        string solver_type_ = "Klu";
        
        vector<Amesos_BaseSolver*> solver_;
        Amesos factory_;
        Teuchos::ParameterList list_;

        int initialize_matrix();
        int initialize_lhs();
        int initialize_rhs();
        int calculate_rhs();
        int initialize_problem();
        int initialize_solver();
        int initialize_transport();
        
        double compute_k1(unsigned cell, unsigned group, unsigned moment);
        double compute_k2(unsigned cell, unsigned group, unsigned moment);
        double compute_k3(unsigned cell, unsigned group, unsigned moment);
    };
}
#endif
