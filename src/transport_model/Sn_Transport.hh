#ifndef Sn_Transport_hh
#define Sn_Transport_hh

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseSolver.h>

#include "Data.hh"
#include "Mesh.hh"
#include "Ordinates.hh"

namespace transport_ns
{
    using std::string;
    using std::vector;
    using std::pow;
    using namespace data_ns;
    using namespace mesh_ns;
    
    class Sn_Transport
    {
    private:

        void calculate_leakage(vector<double> &psi,
                               vector<double> &leakage);

        void calculate_source(vector<double> &q,
                              vector<double> &phi);

        void slab_sweep(vector<double> &psi,
                        vector<double> &q);

        void spherical_sweep(vector<double> &psi,
                             vector<double> &q);

        void sweep_special(vector<double> &psi_half,
                           vector<double> &q);

        void sweep_inward(vector<double> &psi,
                          vector<double> &psi_half,
                          vector<double> &q,
                          vector<double> &psi_boundary_sources);

        void sweep_outward(vector<double> &psi,
                           vector<double> &psi_half,
                           vector<double> &q,
                           vector<double> &psi_boundary_sources);

        void check_convergence(bool &converged,
                               vector<double> &phi,
                               vector<double> &phi_old,
                               vector<double> &error_phi,
                               vector<double> &error_phi_old);

        void epetra_solve(Epetra_SerialDenseMatrix &matrix, Epetra_SerialDenseVector &lhs, Epetra_SerialDenseVector &rhs);        
        
        Data &data_;
        Mesh &mesh_;
        Ordinates &ordinates_;
        
        unsigned max_iterations_ = 2;
        double tolerance_ = 1e-10;
        
        double iterations_;
        vector<double> phi_;

        double get_k_spec(unsigned n1, unsigned n2)
        {
            unsigned o = n2 + mesh_.number_of_nodes() * n1;
            
            switch(o)
            {
            case 0:
                return -1.0 / 2.0;
            case 1:
                return -1.0 / 2.0;
            case 2:
                return 1.0 / 2.0;
            case 3:
                return 1.0 / 2.0;
            default:
                cerr << "no such k_spec" << endl;
                return 0;
            }
        }

        double get_m_spec(unsigned n1, unsigned n2)
        {
            unsigned o = n2 + mesh_.number_of_nodes() * n1;
            
            switch(o)
            {
            case 0:
                return 2.0 / 3.0;
            case 1:
                return 1.0 / 3.0;
            case 2:
                return 1.0 / 3.0;
            case 3:
                return 2.0 / 3.0;
            default:
                cerr << "no such m_spec" << endl;
                return 0;
            }
        }

        double get_k(unsigned i, unsigned n1, unsigned n2)
        {
            unsigned o = n2 + mesh_.number_of_nodes() * n1;
            
            switch(o)
            {
            case 0:
                return - pow(mesh_.cell_center_position(i), 2) / 2.0 + mesh_.cell_center_position(i) * mesh_.cell_length(i) / 6.0 - pow(mesh_.cell_length(i), 2) / 24.0;
            case 1:
                return - pow(mesh_.cell_center_position(i), 2) / 2.0 - mesh_.cell_center_position(i) * mesh_.cell_length(i) / 6.0 - pow(mesh_.cell_length(i), 2) / 24.0;
            case 2:
                return pow(mesh_.cell_center_position(i), 2) / 2.0 - mesh_.cell_center_position(i) * mesh_.cell_length(i) / 6.0 + pow(mesh_.cell_length(i), 2) / 24.0;
            case 3:
                return pow(mesh_.cell_center_position(i), 2) / 2.0 + mesh_.cell_center_position(i) * mesh_.cell_length(i) / 6.0 + pow(mesh_.cell_length(i), 2) / 24.0;
            default:
                cerr << "no such k" << endl;
                return 0;
            }
        }

        double get_l(unsigned i, unsigned n1, unsigned n2)
        {
            unsigned o = n2 + mesh_.number_of_nodes() * n1;
            
            switch(o)
            {
            case 0:
                return mesh_.cell_center_position(i) * 2.0 / 3.0 - mesh_.cell_length(i) / 6.0;
            case 1:
                return mesh_.cell_center_position(i) / 3.0;
            case 2:
                return mesh_.cell_center_position(i) / 3.0;
            case 3:
                return mesh_.cell_center_position(i) * 2.0 / 3.0 + mesh_.cell_length(i) / 6.0;
            default:
                cerr << "no such k" << endl;
                return 0;
            }
        }

        double get_m(unsigned i, unsigned n1, unsigned n2)
        {
            unsigned n = n2 + mesh_.number_of_nodes() * n1;
            
            switch(n)
            {
            case 0:
                return pow(mesh_.cell_center_position(i), 2) * 2.0 / 3.0 - mesh_.cell_center_position(i) * mesh_.cell_length(i) / 3.0 + pow(mesh_.cell_length(i), 2) / 15.0;
            case 1:
                return pow(mesh_.cell_center_position(i), 2) / 3.0 + pow(mesh_.cell_length(i), 2) / 60.0;
            case 2:
                return pow(mesh_.cell_center_position(i), 2) / 3.0 + pow(mesh_.cell_length(i), 2) / 60.0;
            case 3:
                return pow(mesh_.cell_center_position(i), 2) * 2.0 / 3.0 + mesh_.cell_center_position(i) * mesh_.cell_length(i) / 3.0 + pow(mesh_.cell_length(i), 2) / 15.0;
            default:
                cerr << "no such k" << endl;
                return 0;
            }
        }

        double a_pos(unsigned i, unsigned g, unsigned o, unsigned n1, unsigned n2)
        {
            double val = - 2 * ordinates_.ordinates(o) / mesh_.cell_length(i) * get_k(i, n1, n2) + 4 * ordinates_.alpha_half(o) / ordinates_.weights(o) * get_l(i, n1, n2) + data_.sigma_t(i, g) * get_m(i, n1, n2);
            
            if (n1 == 1 && n2 == 1)
            {
                val += 2 * ordinates_.ordinates(o) / mesh_.cell_length(i) * pow(mesh_.cell_edge_position(i, 1), 2);
            }
            
            return val;
        }
        
        double a_neg(unsigned i, unsigned g, unsigned o, unsigned n1, unsigned n2)
        {
            double val = - 2 * ordinates_.ordinates(o) / mesh_.cell_length(i) * get_k(i, n1, n2) + 4 * ordinates_.alpha_half(o) / ordinates_.weights(o) * get_l(i, n1, n2) + data_.sigma_t(i, g) * get_m(i, n1, n2);
            
            if (n1 == 0 && n2 == 0)
            {
                val -= 2 * ordinates_.ordinates(o) / mesh_.cell_length(i) * pow(mesh_.cell_edge_position(i, 0), 2);
            }
            
            return val;
        }
        
        double a_spec(unsigned i, unsigned g, unsigned n1, unsigned n2)
        {
            double val = 2 / mesh_.cell_length(i) * get_k_spec(n1, n2) + data_.sigma_t(i, g) * get_m_spec(n1, n2);
            
            if (n1 == 0 && n2 == 0)
            {
                val += 2 / mesh_.cell_length(i);
            }
            
            return val;
        }
        
    public:
        
        Sn_Transport(Data &data,
                     Mesh &mesh,
                     Ordinates &ordinates);
        
        void solve();
        
        void psi_to_phi(vector<double> &phi,
                        vector<double> &psi);

        void update_psi_half(vector<double> &psi_half,
                             vector<double> &psi,
                             unsigned o);
        
        void print_scalar_flux()
        {
            using namespace std;

            const int w = 8;
            
            cout << "Sn_Transport" << endl;

            cout << left;
            cout << setw(w) <<  "cell" << setw(w) << "node" << setw(w) << "group" << setw(w) << "phi" << endl;
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                for (unsigned i = 0; i < mesh_.number_of_cells(); ++i)
                {
                    for (unsigned n = 0; n < mesh_.number_of_nodes(); ++n)
                    {
                        unsigned k = n + mesh_.number_of_nodes() * (g + data_.number_of_groups() * i);
                        
                        cout << setw(w) << i << setw(w) << n << setw(w) << g << setw(w) << phi_[k] << endl;
                    }
                }
            }

            cout << "it: " << iterations_ << endl << endl;
        }

        void print_angular_flux(vector<double> &psi)
        {
            using namespace std;

            const int w = 10;
            
            cout << "Sn_Transport" << endl;

            cout << left;
            cout << setw(w) <<  "cell" << setw(w) << "node" << setw(w) << "group" << setw(w) << "ordinate" << setw(w) << "phi" << endl;
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                for (unsigned i = 0; i < mesh_.number_of_cells(); ++i)
                {
                    for (unsigned n = 0; n < mesh_.number_of_nodes(); ++n)
                    {
                        for (unsigned o = 0; o < ordinates_.number_of_ordinates(); ++o)
                        {
                            unsigned k = n + mesh_.number_of_nodes() * (g + data_.number_of_groups() * (o + ordinates_.number_of_ordinates() * i));
                            
                            cout << setw(w) << i << setw(w) << n << setw(w) << g << setw(w) << o << setw(w) << psi[k] << endl;
                        }
                    }
                }
            }

            cout << "it: " << iterations_ << endl << endl;
        }
    };
}
#endif
