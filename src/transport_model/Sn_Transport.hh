#ifndef Lumped_Linear_Discontinuous_hh
#define Lumped_Linear_Discontinuous_hh

#include <iostream>
#include <string>
#include <vector>

#include "Data.hh"
#include "Mesh.hh"
#include "Ordinates.hh"

namespace transport_ns
{
    using std::string;
    using std::vector;
    using namespace data_ns;
    using namespace mesh_ns;
    
    class Sn_Transport
    {
    private:

        // vector operators
    
        void calculate_leakage(vector<double> &psi,
                               vector<double> &leakage);

        void calculate_source(vector<double> &q,
                              vector<double> &phi);

        void lumped_linear_discontinuous_sweep(vector<double> &psi,
                                               vector<double> &q);
    
        void check_convergence(bool &converged,
                               vector<double> &phi,
                               vector<double> &phi_old,
                               vector<double> &error_phi,
                               vector<double> &error_phi_old);
    
    
        // data

        Data &data_;
        Mesh &mesh_;
        Ordinates &ordinates_;
        
        unsigned max_iterations_ = 1000;
        double tolerance_ = 1e-8;
        
        string problem_type_;

        double iterations_;
        vector<double> phi_;

    public:

        Sn_Transport(Data &data,
                     Mesh &mesh,
                     Ordinates &ordinates,
                     string problem_type = "forward");

        void solve(vector<double> &psi,
                   vector<double> &leakage,
                   unsigned &iterations);

        void psi_to_phi(vector<double> &phi,
                        vector<double> &psi);
    
        void print_scalar_flux()
        {
            std::cout << "Sn_Transport" << std::endl;
            
            for (unsigned i = 0; i < mesh_.number_of_cells(); ++i)
            {
                for (unsigned g = 0; g < mesh_.number_of_cells(); ++g)
                {
                    for (unsigned n = 0; n < mesh_.number_of_nodes(); ++n)
                    {
                        unsigned k = n + mesh_.number_of_nodes() * (g + data_.number_of_groups() * i);
                        
                        std::cout << phi_[k] << std::endl;
                }
            }
        }
    };
}
#endif
