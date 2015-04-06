#include "Sn_Transport.hh"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <mpi.h>

#include "Data.hh"
#include "Mesh.hh"
#include "Ordinates.hh"

namespace transport_ns
{
    using namespace std;
    using namespace data_ns;
    using namespace mesh_ns;
    
    Sn_Transport::
    Sn_Transport(Data &data,
                 Mesh &mesh,
                 Ordinates &ordinates,
                 string problem_type)
    :
        data_(data),
        mesh_(mesh),
        ordinates_(ordinates),
        problem_type_(problem_type)
    {
        data_.check();
        mesh_.check();
    }
    
    void Sn_Transport::
    calculate_source(vector<double> &q,
                     vector<double> &phi)
    {
        // calculate the source
        for (unsigned i = 0; i < mesh_.number_of_cells(); ++i)
        {
            for (unsigned n = 0; n < mesh_.number_of_nodes(); ++n)
            {
                for (unsigned g = 0; g < data_.number_of_groups(); ++g)
                {
                    double sum = 0;
                    for (unsigned g2 = 0; g2 < data_.number_of_groups(); ++g2)
                    {
                        unsigned k3 = n + mesh_.number_of_nodes() * (g2 + data_.number_of_groups() * i);
                        
                        sum += (data_.chi(i, g) * data_.nu_sigma_f(i, g2) + data_.sigma_s(i, g2, g)) * phi[k3];
                    }
                    unsigned k1 = n + mesh_.number_of_nodes() * (g + data_.number_of_groups() * i);
                    
                    q[k1] = (data_.internal_source(i, g) + sum) / 2;
                }
            }
        }
    }

    void Sn_Transport::
    lumped_linear_discontinuous_sweep(vector<double> &psi,
                                      vector<double> &q)
    {
        vector<double> psi_boundary_sources(data_.number_of_groups() * ordinates_.number_of_ordinates(), 0);
    
        // boundary condition, x=0
        if (boundary_conditions[0] == "reflected")
        {
            for (unsigned o = 0; o < ordinates_.number_of_ordinates() / 2; ++o)
            {
                for (unsigned g = 0; g < data_.number_of_groups(); ++g)
                {
                    unsigned k1 = g + data_.number_of_groups() * o;
                    unsigned k2 = 0 + mesh_.number_of_nodes() * (g + data_.number_of_groups() * ((ordinates_.number_of_ordinates() - 1 - o) + ordinates_.number_of_ordinates() * 0));
                    
                    psi_boundary_sources[k1] = psi[k2];
                }
            }
        }
        else if (boundary_conditions[0] == "vacuum")
        {
            for (unsigned o = 0; o < ordinates_.number_of_ordinates() / 2; ++o)
            {
                for (unsigned g = 0; g < data_.number_of_groups(); ++g)
                {
                    unsigned k1 = g + data_.number_of_groups() * o;
                                        
                    psi_boundary_sources[k1] = 0;
                }
            }
        }
        
        // sweep first cell
        for (unsigned o = 0; o < ordinates_.number_of_ordinates() / 2; ++o)
        {
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                unsigned i = 0;
                unsigned k3 = mesh_.number_of_nodes() * (g + data_.number_of_groups() * i);
                unsigned k4 = mesh_.number_of_nodes() * (g + data_.number_of_groups() * (o + ordinates_.number_of_ordinates() * i));
				
                double a1 = ordinates_.ordinates(o) / 2 + data_.sigma_t(i, g) * mesh_.cell_length(i) / 2;
                double a2 = ordinates_.ordinates(o) / 2;
                double a3 = -ordinates_.ordinates(o) / 2;
                double a4 = ordinates_.ordinates(o) / 2 + data_.sigma_t(i, g) * mesh_.cell_length(i) / 2;
                double s1 = mesh_.cell_length(i) * (q[k3] / 2) + ordinates_.ordinates(o) * (data_.boundary_source(0, g, o) + psi_boundary_sources[k2]);
                double s2 = mesh_.cell_length(i) * (q[k3+1] / 2);
				
                psi[k4] = (a4 * s1 - a2 * s2) / (a1 * a4 - a2 * a3);
                psi[k4+1] = (a3 * s1 - a1 * s2) / (a2 * a3 - a1 * a4);
            }
        }

        // sweep right over cells
        for (unsigned i = 1; i < mesh_.number_of_cells(); ++i)
        {
            for (unsigned o = 0; o < ordinates_.number_of_ordinates() / 2; ++o)
            {
                for (unsigned g = 0; g < data_.number_of_groups(); ++g)
                {
                    unsigned k2 = 1 + mesh_.number_of_nodes() * (g + data_.number_of_groups() * (o + ordinates_.number_of_ordinates() * (i-1)));
                    unsigned k3 = mesh_.number_of_nodes() * (g + data_.number_of_groups() * i);
                    unsigned k4 = mesh_.number_of_nodes() * (g + data_.number_of_groups() * (o + ordinates_.number_of_ordinates() * i));
				
                    double a1 = ordinates_.ordinates(o) / 2 + data_.sigma_t(i, g) * mesh_.cell_length(i) / 2;
                    double a2 = ordinates_.ordinates(o) / 2;
                    double a3 = -ordinates_.ordinates(o) / 2;
                    double a4 = ordinates_.ordinates(o) / 2 + data_.sigma_t(i, g) * mesh_.cell_length(i) / 2;
                    double s1 = mesh_.cell_length(i) * (q[k3] / 2) + ordinates_.ordinates(o) * psi[k2];
                    double s2 = mesh_.cell_length(i) * (q[k3+1] / 2);
				
                    psi[k4] = (a4 * s1 - a2 * s2) / (a1 * a4 - a2 * a3);
                    psi[k4+1] = (a3 * s1 - a1 * s2) / (a2 * a3 - a1 * a4);
                }
            }
        }

        // boundary condition, x=X
        if (boundary_conditions[1] == "reflected")
        {
            for (unsigned o = ordinates_.number_of_ordinates() / 2; o < ordinates_.number_of_ordinates(); ++o)
            {
                for (unsigned g = 0; g < data_.number_of_groups(); ++g)
                {
                    unsigned k1 = g + data_.number_of_groups() * o;
                    unsigned k2 = 1 + mesh_.number_of_nodes() * (g + data_.number_of_groups() * ((ordinates_.number_of_ordinates() - 1 - o) + ordinates_.number_of_ordinates() * (mesh_.number_of_cells() - 1)));
                    
                    psi_boundary_sources[k1] = psi[k2];
                }
            }
        }
        else if (boundary_conditions[1] == "vacuum")
        {
            for (unsigned o = ordinates_.number_of_ordinates() / 2; o < ordinates_.number_of_ordinates(); ++o)
            {
                for (unsigned g = 0; g < data_.number_of_groups(); ++g)
                {
                    unsigned k1 = g + data_.number_of_groups() * o;
                                        
                    psi_boundary_sources[k1] = 0;
                }
            }
        }


        // sweep final cell
        for (unsigned o = ordinates_.number_of_ordinates() / 2; o < ordinates_.number_of_ordinates(); ++o)
        {
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                unsigned i = mesh_.number_of_cells() - 1;
                unsigned k3 = mesh_.number_of_nodes() * (g + data_.number_of_groups() * i);
                unsigned k4 = mesh_.number_of_nodes() * (g + data_.number_of_groups() * (o + ordinates_.number_of_ordinates() * i));
				
                double a1 = -ordinates_.ordinates(o) / 2 + data_.sigma_t(i, g) * mesh_.cell_length(i) / 2;
                double a2 = ordinates_.ordinates(o) / 2;
                double a3 = -ordinates_.ordinates(o) / 2;
                double a4 = -ordinates_.ordinates(o) / 2 + data_.sigma_t(i, g) * mesh_.cell_length(i) / 2;
                double s1 = mesh_.cell_length(i) * (q[k3] / 2);
                double s2 = mesh_.cell_length(i) * (q[k3+1] / 2) - ordinates_.ordinates(o) * (data_.boundary_source(1, g, o) + psi_boundary_sources[k2]);
				
                psi[k4] = (a4 * s1 - a2 * s2) / (a1 * a4 - a2 * a3);
                psi[k4+1] = (a3 * s1 - a1 * s2) / (a2 * a3 - a1 * a4);		
            }
        }

        // sweep left over cells
        for (int i = mesh_.number_of_cells() - 2; i >= 0; --i)
        {
            for (unsigned o = ordinates_.number_of_ordinates() / 2; o < ordinates_.number_of_ordinates(); ++o)
            {
                for (unsigned g = 0; g < data_.number_of_groups(); ++g)
                {
                    unsigned k2 = mesh_.number_of_nodes() * (g + data_.number_of_groups() * (o + ordinates_.number_of_ordinates() * (i+1))); // n=1
                    unsigned k3 = mesh_.number_of_nodes() * (g + data_.number_of_groups() * i);
                    unsigned k4 = mesh_.number_of_nodes() * (g + data_.number_of_groups() * (o + ordinates_.number_of_ordinates() * i));
				
                    double a1 = -ordinates_.ordinates(o) / 2 + data_.sigma_t(i, g) * mesh_.cell_length(i) / 2;
                    double a2 = ordinates_.ordinates(o) / 2;
                    double a3 = -ordinates_.ordinates(o) / 2;
                    double a4 = -ordinates_.ordinates(o) / 2 + data_.sigma_t(i, g) * mesh_.cell_length(i) / 2;
                    double s1 = mesh_.cell_length(i) * (q[k3] / 2);
                    double s2 = mesh_.cell_length(i) * (q[k3+1] / 2) - ordinates_.ordinates(o) * psi[k2];
                    
                    psi[k4] = (a4 * s1 - a2 * s2) / (a1 * a4 - a2 * a3);
                    psi[k4+1] = (a3 * s1 - a1 * s2) / (a2 * a3 - a1 * a4);		
                }
            }
        }
    }

    void Sn_Transport::
    psi_to_phi(vector<double> &phi,
               vector<double> &psi)
    {
        for (unsigned i = 0; i < mesh_.number_of_cells(); ++i)
        {
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                for (unsigned n = 0; n < mesh_.number_of_nodes(); ++n)
                {
                    double sum = 0;
		
                    for (unsigned o = 0; o < ordinates_.number_of_ordinates(); ++o)
                    {
                        unsigned k1 = n + mesh_.number_of_nodes() * (g + data_.number_of_groups() * (o + ordinates_.number_of_ordinates() * i));
			
                        sum += ordinates_.weights(o) * psi[k1];
                    }
                    unsigned k1 = n + mesh_.number_of_nodes() * (g + data_.number_of_groups() * i);
                    
                    phi[k1] = sum;
                }
            }
        }
    }

    void Sn_Transport::
    check_convergence(bool &converged,
                      vector<double> &phi,
                      vector<double> &phi_old,
                      vector<double> &error_phi,
                      vector<double> &error_phi_old)
    {
        error_phi_old = error_phi;
    
        converged = true;
        
        for (unsigned i = 0; i < mesh_.number_of_cells(); ++i)
        {
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                for (unsigned n = 0; n < mesh_.number_of_nodes(); ++n)
                {
                    unsigned k1 = n + mesh_.number_of_nodes() * (g + data_.number_of_groups() * i);
                    
                    error_phi[k1] = abs(phi[k1] - phi_old[k1]) / phi[k1];
                    
                    if (error_phi[k1] > tolerance_ * (1 - error_phi[k1] / error_phi_old[k1]))
                    {
                        converged = false;
                    }
                }
            }
        }
    }

    void Sn_Transport::
    solve()
    {
        using namespace std;
        
        // temporary variables
        vector<double> psi(mesh_.number_of_cells() * data_.number_of_groups() * mesh_.number_of_nodes() * ordinates_.number_of_ordinates());
        vector<double> phi(mesh_.number_of_cells() * data_.number_of_groups() * mesh_.number_of_nodes(), 0); // average scalar flux
        vector<double> phi_old(mesh_.number_of_cells() * data_.number_of_groups() * mesh_.number_of_nodes(), 0); // average scalar flux from previous iteration
        vector<double> q(mesh_.number_of_cells() * data_.number_of_groups() * mesh_.number_of_nodes(), 0); // total source, including fission and scattering
        vector<double> error_phi(mesh_.number_of_cells() * data_.number_of_groups() * mesh_.number_of_nodes(), 1);
        vector<double> error_phi_old(mesh_.number_of_cells() * data_.number_of_groups() * mesh_.number_of_nodes(), 1);
        bool converged = false;
    
        // begin iterations
        for (unsigned it = 0; it < max_iterations_; ++it)
        {
            calculate_source(q,
                             phi);
        
            lumped_linear_discontinuous_sweep(psi,
                                              q);
            
            phi_old = phi;
            psi_to_phi(phi,
                       psi);
            check_convergence(converged,
                              phi,
                              phi_old,
                              error_phi,
                              error_phi_old);
            
            if (converged)
            {
                iterations_ = it + 1;
                phi_ = phi;
                break;
            }
            else if (it==max_iterations_ - 1)
            {
                iterations = it + 1;
            }
        }
        // calculate_leakage(psi,
        //                   leakage);
    }

    void Sn_Transport::
    calculate_leakage(vector<double> &psi,
                      vector<double> &leakage)
    {
        unsigned il = 0;
        unsigned ir = mesh_.number_of_cells() - 1;

        leakage.assign(leakage.size(), 0.0);
        
        for (unsigned op = 0; op < ordinates_.number_of_ordinates() / 2; ++op)
        {
            unsigned on = ordinates_.number_of_ordinates() - op - 1;
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                unsigned kl = 0 + mesh_.number_of_nodes() * (g + data_.number_of_groups() * (on + ordinates_.number_of_ordinates() * il));
                unsigned kr = 1 + mesh_.number_of_nodes() * (g + data_.number_of_groups() * (op + ordinates_.number_of_ordinates() * ir));
                
                leakage[0] += ordinates_.ordinates(on) * ordinates_.weights(on) * psi[kl];
                leakage[1] += ordinates_.ordinates(op) * ordinates_.weights(op) * psi[kr];
            }
        }
    }
}
