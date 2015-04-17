#ifndef FEM_Transport_hh
#define FEM_Transport_hh

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "Data.hh"
#include "Mesh.hh"

namespace transport_ns
{
    using std::string;
    using std::vector;
    using std::pow;
    using namespace data_ns;
    using namespace mesh_ns;
    
    class FEM_Transport
    {
    private:

        void calculate_source();

        void spherical_sweep();

        void sweep_inward();

        void sweep_outward();
        
        void check_convergence();

        void analytic_solve(vector<double> &matrix, vector<double> &lhs, vector<double> &rhs);

        void epetra_solve(vector<double> &matrix, vector<double> &lhs, vector<double> &rhs, int size);
        
        Data &data_;
        Mesh &mesh_;

        bool converged_;
        
        unsigned number_of_angular_nodes_ = 2;
        unsigned total_nodes_ = 4;
        unsigned max_iterations_;
        double tolerance_;

        vector<double> xi_coord = {-1.0, -1.0, 1.0, 1.0};
        vector<double> eta_coord = {-1.0, 1.0, -1.0, 1.0};
        vector<double> mp_coord = {-1.0, 1.0};
        
        unsigned iterations_;
        vector<double> phi_;
        vector<double> phi_old_;
        vector<double> error_phi_;
        vector<double> error_phi_old_;
        vector<double> psi_;
        vector<double> q_;
        vector<double> psi_boundary_sources_;
        double k_eigenvalue_;
        //double k_eigenvalue_old_;

        vector<double> matrix_;
        vector<double> rhs_;
        vector<double> lhs_;
        
        double get_k_r(unsigned i, unsigned n1, unsigned n2, unsigned mp);

        double get_k_l(unsigned i, unsigned n1, unsigned n2, unsigned mp);

        double get_j(unsigned i, unsigned n1, unsigned n2, unsigned mp);
        
        double get_l(unsigned i, unsigned n1, unsigned n2, unsigned mp);

        double get_m(unsigned i, unsigned n1, unsigned n2);

        inline double &psi(unsigned i, unsigned g, unsigned n, unsigned o, unsigned mp)
        {
            return psi_[n + mesh_.number_of_nodes() * (g + data_.number_of_groups() * (o + number_of_angular_nodes_ * (mp + 2 * i)))];
        }
        
        inline double &q(unsigned i, unsigned g, unsigned n)
        {
            return q_[n + mesh_.number_of_nodes() * (g + data_.number_of_groups() * i)];
        }

        inline double &phi(unsigned i, unsigned g, unsigned n)
        {
            return phi_[n + mesh_.number_of_nodes() * (g + data_.number_of_groups() * i)];
        }

        inline double &phi_old(unsigned i, unsigned g, unsigned n)
        {
            return phi_old_[n + mesh_.number_of_nodes() * (g + data_.number_of_groups() * i)];
        }

        inline double &error_phi(unsigned i, unsigned g, unsigned n)
        {
            return error_phi_[n + mesh_.number_of_nodes() * (g + data_.number_of_groups() * i)];
        }
        
        inline double &error_phi_old(unsigned i, unsigned g, unsigned n)
        {
            return error_phi_old_[n + mesh_.number_of_nodes() * (g + data_.number_of_groups() * i)];
        }
        
        inline double &psi_boundary_sources(unsigned g, unsigned o, unsigned mp)
        {
            return psi_boundary_sources_[g + data_.number_of_groups() * (o + number_of_angular_nodes_ * mp)];
        }

        inline double &matrix(unsigned al, unsigned be)
        {
            return matrix_[be + number_of_angular_nodes_ * mesh_.number_of_nodes() * al];
        }

        inline double &lhs(unsigned n, unsigned o)
        {
            return lhs(o + number_of_angular_nodes_ * n);
        }

        inline double &lhs(unsigned be)
        {
            return lhs_[be];
        }
        
        inline double &rhs(unsigned n, unsigned o)
        {
            return rhs(n + mesh_.number_of_nodes() * o);
        }

        inline double &rhs(unsigned be)
        {
            return rhs_[be];
        }

        
    public:
        
        FEM_Transport(Data &data,
                      Mesh &mesh,
                      unsigned max_iterations,
                      double tolerance);
        
        void solve();
        
        void solve_eigenvalue();
        
        void psi_to_phi();

        void print_eigenvalue()
        {
            cout << "FEM_Transport" << endl;
            
            cout << "k_eigenvalue = " << setprecision(15) << k_eigenvalue_ << endl;

            cout << "iterations: " << iterations_ << endl << endl;
        }

        void calculate_k();

        void normalize_phi();
        
        void print_scalar_flux()
        {
            using namespace std;

            const int w = 8;
            
            cout << "FEM_Transport" << endl;

            cout << left;
            cout << setw(w) <<  "cell" << setw(w) << "node" << setw(w) << "group" << setw(w) << "phi" << endl;
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                for (unsigned i = 0; i < mesh_.number_of_cells(); ++i)
                {
                    for (unsigned n = 0; n < mesh_.number_of_nodes(); ++n)
                    {
                        cout << setw(w) << i << setw(w) << n << setw(w) << g << setw(w) << phi(i, g, n) << endl;
                    }
                }
            }

            cout << "iterations: " << iterations_ << endl << endl;
        }

        void print_angular_flux()
        {
            using namespace std;

            const int w = 10;
            
            cout << "FEM_Transport" << endl;

            cout << left;
            cout << setw(w) <<  "cell" << setw(w) << "node" << setw(w) << "group" << setw(w) << "ordinate" << setw(w) << "neg" << setw(w) << "phi" << endl;
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                for (unsigned i = 0; i < mesh_.number_of_cells(); ++i)
                {
                    for (unsigned n = 0; n < mesh_.number_of_nodes(); ++n)
                    {
                        for (unsigned o = 0; o < number_of_angular_nodes_; ++o)
                        {
                            for (unsigned mp = 0; mp < 2; ++mp)
                            {
                                cout << setw(w) << i << setw(w) << n << setw(w) << g << setw(w) << o << setw(w) << mp << setw(w) << psi(i, g, n, o, mp) << endl;
                            }
                        }
                    }
                }
            }

            cout << "it: " << iterations_ << endl << endl;
        }
    };
}
#endif
