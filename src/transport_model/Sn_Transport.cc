#include "Sn_Transport.hh"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <mpi.h>

#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseSolver.h>

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
                 Ordinates &ordinates)
    :
        data_(data),
        mesh_(mesh),
        ordinates_(ordinates)
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
    slab_sweep(vector<double> &psi,
               vector<double> &q)
    {
        vector<double> psi_boundary_sources(data_.number_of_groups() * ordinates_.number_of_ordinates(), 0);
    
        // boundary condition, x=0
        if (data_.boundary_condition(0).compare("reflected") == 0)
        {
            for (unsigned o = ordinates_.number_of_ordinates() / 2; o < ordinates_.number_of_ordinates(); ++o)
            {
                for (unsigned g = 0; g < data_.number_of_groups(); ++g)
                {
                    unsigned k1 = g + data_.number_of_groups() * o;
                    unsigned k2 = 0 + mesh_.number_of_nodes() * (g + data_.number_of_groups() * ((ordinates_.number_of_ordinates() - 1 - o) + ordinates_.number_of_ordinates() * 0));
                    
                    psi_boundary_sources[k1] = psi[k2];
                }
            }
        }
        else if (data_.boundary_condition(0).compare("vacuum") == 0)
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
        else
        {
            cerr << "boundary condition \"" << data_.boundary_condition(0) << "\" not available" << endl;
        }
        
        // sweep first cell
        for (unsigned o = ordinates_.number_of_ordinates() / 2; o < ordinates_.number_of_ordinates(); ++o)
        {
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                unsigned i = 0;
                unsigned k2 = g + data_.number_of_groups() * o;
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
            for (unsigned o = ordinates_.number_of_ordinates() / 2; o < ordinates_.number_of_ordinates(); ++o)
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
        if (data_.boundary_condition(1).compare("reflected") == 0)
        {
            for (unsigned o = 0; o < ordinates_.number_of_ordinates() / 2; ++o)
            {
                for (unsigned g = 0; g < data_.number_of_groups(); ++g)
                {
                    unsigned k1 = g + data_.number_of_groups() * o;
                    unsigned k2 = 1 + mesh_.number_of_nodes() * (g + data_.number_of_groups() * ((ordinates_.number_of_ordinates() - 1 - o) + ordinates_.number_of_ordinates() * (mesh_.number_of_cells() - 1)));
                    
                    psi_boundary_sources[k1] = psi[k2];
                }
            }
        }
        else if (data_.boundary_condition(1).compare("vacuum") == 0)
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
        else
        {
            cerr << "boundary condition \"" << data_.boundary_condition(1) << "\" not available" << endl;
        }
        
        // sweep final cell
        for (unsigned o = 0; o < ordinates_.number_of_ordinates() / 2; ++o)
        {
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                unsigned i = mesh_.number_of_cells() - 1;
                unsigned k2 = g + data_.number_of_groups() * o;
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
            for (unsigned o = 0; o < ordinates_.number_of_ordinates() / 2; ++o)
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
    spherical_sweep(vector<double> &psi,
                    vector<double> &q)
    {
        vector<double> psi_boundary_sources(data_.number_of_groups() * ordinates_.number_of_ordinates(), 0);
        vector<double> psi_half(mesh_.number_of_nodes() * data_.number_of_groups() * mesh_.number_of_cells());
        
        // boundary condition, r=R
        if (data_.boundary_condition(1).compare("reflected") == 0)
        {
            for (unsigned o = 0; o < ordinates_.number_of_ordinates() / 2; ++o)
            {
                for (unsigned g = 0; g < data_.number_of_groups(); ++g)
                {
                    unsigned k1 = g + data_.number_of_groups() * o;
                    unsigned k2 = 1 + mesh_.number_of_nodes() * (g + data_.number_of_groups() * ((ordinates_.number_of_ordinates() - 1 - o) + ordinates_.number_of_ordinates() * (mesh_.number_of_cells() - 1)));
                    
                    psi_boundary_sources[k1] = psi[k2] + data_.boundary_source(1, g, o);
                }
            }
        }
        else if (data_.boundary_condition(1).compare("vacuum") == 0)
        {
            for (unsigned o = 0; o < ordinates_.number_of_ordinates() / 2; ++o)
            {
                for (unsigned g = 0; g < data_.number_of_groups(); ++g)
                {
                    unsigned k = g + data_.number_of_groups() * o;
                    
                    psi_boundary_sources[k] = data_.boundary_source(1, g, o);
                }
            }

        }
        else
        {
            cerr << "boundary condition \"" << data_.boundary_condition(1) << "\" not available" << endl;
        }

        // sweep inward
        
        sweep_special(psi_half, q);
        
        sweep_inward(psi, psi_half, psi_boundary_sources, q);
        
        // boundary condition, r=0
        if (data_.boundary_condition(0).compare("reflected") == 0)
        {
            for (unsigned o = ordinates_.number_of_ordinates() / 2; o < ordinates_.number_of_ordinates(); ++o)
            {
                for (unsigned g = 0; g < data_.number_of_groups(); ++g)
                {
                    unsigned k1 = g + data_.number_of_groups() * o;
                    unsigned k2 = 0 + mesh_.number_of_nodes() * (g + data_.number_of_groups() * ((ordinates_.number_of_ordinates() - 1 - o) + ordinates_.number_of_ordinates() * 0));
                    
                    psi_boundary_sources[k1] = psi[k2];
                }
            }
        }
        else
        {
            cerr << "boundary condition \"" << data_.boundary_condition(0) << "\" not available" << endl;
        }
        
        // sweep outward
        
        sweep_outward(psi, psi_half, psi_boundary_sources, q);
    }

    void Sn_Transport::
    sweep_special(vector<double> &psi_half, vector<double> &q)
    {
        for (int i = mesh_.number_of_cells() - 1; i >= 0; --i)
        {
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                Epetra_SerialDenseMatrix matrix (mesh_.number_of_nodes(), mesh_.number_of_nodes());
                Epetra_SerialDenseVector lhs (mesh_.number_of_nodes());
                Epetra_SerialDenseVector rhs (mesh_.number_of_nodes());
        
                for (unsigned n1 = 0; n1 < mesh_.number_of_nodes(); ++n1)
                {
                    double sum = 0;
                    
                    for (unsigned n2 = 0; n2 < mesh_.number_of_nodes(); ++n2)
                    {
                        matrix(n1, n2) = a_spec(i, g, n1, n2);
                        
                        unsigned k = mesh_.number_of_nodes() * (g + data_.number_of_groups() * i);
                        
                        sum += get_m_spec(n1, n2) * q[k] / 2;
                    }
                    
                    rhs(n1) = sum;
                }
                
                if (i < mesh_.number_of_cells() - 1)
                {
                    unsigned k = mesh_.number_of_nodes() * (g + data_.number_of_groups() * (i + 1));
                    
                    rhs(1) += 2.0 / mesh_.cell_length(i) * psi_half[k];
                }
                
                epetra_solve(matrix, lhs, rhs);
                
                for (unsigned n = 0; n < mesh_.number_of_nodes(); ++n)
                {
                    unsigned k = n + mesh_.number_of_nodes() * (g + data_.number_of_groups() * i);
                    
                    psi_half[k] = lhs[n];
                }

                // cout << "************************************************" << endl;
                // cout << "special, i:" << i << " g:" << g << endl;
                // cout << lhs << endl;
            }
        }
    }
    
    void Sn_Transport::
    sweep_inward(vector<double> &psi, vector<double> &psi_half, vector<double> &psi_boundary_sources, vector<double> &q)
    {
        for (unsigned o = 0; o < ordinates_.number_of_ordinates() / 2; ++o)
        {
            for (int i = mesh_.number_of_cells() - 1; i >= 0; --i)
            { 
                for (unsigned g = 0; g < data_.number_of_groups(); ++g)
                {
                    Epetra_SerialDenseMatrix matrix (mesh_.number_of_nodes(), mesh_.number_of_nodes());
                    Epetra_SerialDenseVector lhs (mesh_.number_of_nodes());
                    Epetra_SerialDenseVector rhs (mesh_.number_of_nodes());
        
                    for (unsigned n1 = 0; n1 < mesh_.number_of_nodes(); ++n1)
                    {
                        double sum1 = 0;
                        double sum2 = 0;
                        
                        for (unsigned n2 = 0; n2 < mesh_.number_of_nodes(); ++n2)
                        {
                            matrix(n1, n2) = a_neg(i, g, o, n1, n2);

                            unsigned k1 = n2 + mesh_.number_of_nodes() * (g + data_.number_of_groups() * i);
                            unsigned k2 = n2 + mesh_.number_of_nodes() * (g + data_.number_of_groups() * i);

                            sum1 += get_l(i, n1, n2) * psi_half[k1];
                            sum2 += get_m(i, n1, n2) * q[k2] / 2; // isotropic internal source
                        }

                        rhs(n1) = 2 * ordinates_.alpha(o) / ordinates_.weights(o) * sum1 + sum2;
                    }
                    
                    if (i < mesh_.number_of_cells() - 1)
                    {
                        unsigned k = mesh_.number_of_nodes() * (g + data_.number_of_groups() * (i + 1));
                        
                        rhs(1) -= 2 * ordinates_.ordinates(o) / mesh_.cell_length(i) * pow(mesh_.cell_edge_position(i, 1), 2) * psi[k];
                    }
                    else
                    {
                        unsigned k = g + data_.number_of_groups() * o;
                        
                        rhs(1) -= 2 * ordinates_.ordinates(o) / mesh_.cell_length(i) * pow(mesh_.cell_edge_position(i, 1), 2) * psi_boundary_sources[k];
                    }
        
                    epetra_solve(matrix, lhs, rhs);

                    for (unsigned n = 0; n < mesh_.number_of_nodes(); ++n)
                    {
                        unsigned k = n + mesh_.number_of_nodes() * (g + data_.number_of_groups() * (o + ordinates_.number_of_ordinates() * i));
                        
                        psi[k] = lhs[n];
                    }

                    // if (i==1 && o==0)
                    // {
                    //     cout << "************************************************" << endl;
                    //     cout << "inward, i:" << i << " o:" << o << " g:" << g << endl << endl;

                    // cout << "mu_m=" << ordinates_.ordinates(o) << " dri=" << mesh_.cell_length(i) << " aph=" << ordinates_.alpha_half(o) << " a=" << ordinates_.alpha(o) << " r_imh=" << mesh_.cell_edge_position(i, 0) << " w_m=" << ordinates_.weights(o) << " sigma_t=" << data_.sigma_t(i, g) << endl;
                    
                        // cout << matrix << endl;
                        // cout << lhs << endl;
                        // cout << rhs << endl;
                    // }
                }
            }
            
            update_psi_half(psi_half, psi, o);
        }
     }

    void Sn_Transport::
    sweep_outward(vector<double> &psi, vector<double> &psi_half, vector<double> &psi_boundary_sources, vector<double> &q)
    {
        for (unsigned o = ordinates_.number_of_ordinates() / 2; o < ordinates_.number_of_ordinates(); ++o)
        {
            for (unsigned i = 0; i < mesh_.number_of_cells(); ++i)
            { 
                for (unsigned g = 0; g < data_.number_of_groups(); ++g)
                {
                    Epetra_SerialDenseMatrix matrix (mesh_.number_of_nodes(), mesh_.number_of_nodes());
                    Epetra_SerialDenseVector lhs (mesh_.number_of_nodes());
                    Epetra_SerialDenseVector rhs (mesh_.number_of_nodes());
        
                    for (unsigned n1 = 0; n1 < mesh_.number_of_nodes(); ++n1)
                    {
                        double sum1 = 0;
                        double sum2 = 0;
                        
                        for (unsigned n2 = 0; n2 < mesh_.number_of_nodes(); ++n2)
                        {
                            matrix(n1, n2) = a_pos(i, g, o, n1, n2);
                            
                            unsigned k1 = n2 + mesh_.number_of_nodes() * (g + data_.number_of_groups() * i);
                            unsigned k2 = n2 + mesh_.number_of_nodes() * (g + data_.number_of_groups() * i);
                            if (i == 0)
                            {
                                cout << "n1=" << n1 << " n2=" << n2 << " a_pos=" << a_pos(i, g, o, n1, n2) << endl;
                            }
                            sum1 += get_l(i, n1, n2) * psi_half[k1];
                            sum2 += get_m(i, n1, n2) * q[k2] / 2; // isotropic internal source
                        }
                        
                        rhs(n1) = 2 * ordinates_.alpha(o) / ordinates_.weights(o) * sum1 + sum2;
                    }
                    
                    if (i > 0)
                    {
                        unsigned k = mesh_.number_of_nodes() * (g + data_.number_of_groups() * (i - 1));
                        
                        rhs(0) += 2 * ordinates_.ordinates(o) / mesh_.cell_length(i) * pow(mesh_.cell_edge_position(i, 0), 2) * psi[k];
                    }
                    else
                    {
                        unsigned k = g + data_.number_of_groups() * o;
                        
                        rhs(0) += 2 * ordinates_.ordinates(o) / mesh_.cell_length(i) * pow(mesh_.cell_edge_position(i, 0), 2) * psi_boundary_sources[k];
                    }
                    
                    epetra_solve(matrix, lhs, rhs);
                    
                    for (unsigned n = 0; n < mesh_.number_of_nodes(); ++n)
                    {
                        unsigned k = n + mesh_.number_of_nodes() * (g + data_.number_of_groups() * (o + ordinates_.number_of_ordinates() * i));
                        
                        psi[k] = lhs[n];
                    }
                    
                    if (i==0)
                    {
                        cout << "************************************************" << endl;
                        cout << "outward, i:" << i << " o:" << o << " g:" << g << endl << endl;

                        cout << "mu_m=" << ordinates_.ordinates(o) << " dri=" << mesh_.cell_length(i) << " aph=" << ordinates_.alpha_half(o) << " a=" << ordinates_.alpha(o) << " r_imh=" << mesh_.cell_edge_position(i, 0) << " w_m=" << ordinates_.weights(o) << " sigma_t=" << data_.sigma_t(i, g) << endl;
                        
                        cout << matrix << endl;
                        cout << lhs << endl;
                        cout << rhs << endl;
                    }
                }
            }
            
            update_psi_half(psi_half, psi, o);
        }
    }
    
    void Sn_Transport::
    update_psi_half(vector<double> &psi_half,
                    vector<double> &psi,
                    unsigned o)
    {
        for (unsigned i = 0; i < mesh_.number_of_cells(); ++i)
        {
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                for (unsigned n = 0; n < mesh_.number_of_nodes(); ++n)
                {
                    unsigned k1 = n + mesh_.number_of_nodes() * (g + data_.number_of_groups() * i);
                    unsigned k2 = n + mesh_.number_of_nodes() * (g + data_.number_of_groups() * (o + ordinates_.number_of_ordinates() * i));
                    
                    psi_half[k1] = 2 * psi[k2] - psi_half[k1];
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
                    
                    error_phi[k1] = abs((phi[k1] - phi_old[k1]) / phi[k1]);
                    
                    if (error_phi[k1] > tolerance_ * (1 - error_phi[k1] / error_phi_old[k1]))
                    {
                        converged = false;
                    }
                }
            }
        }
    }

    void Sn_Transport::
    epetra_solve(Epetra_SerialDenseMatrix &matrix, Epetra_SerialDenseVector &lhs, Epetra_SerialDenseVector &rhs)
    {
        Epetra_SerialDenseSolver solver;
        
        solver.SetMatrix(matrix);
        solver.SetVectors(lhs, rhs);
        solver.Solve();
    }
    
    void Sn_Transport::
    solve()
    {
        using namespace std;

        if (mesh_.geometry("slab"))
        {
            // temporary variables
            vector<double> psi(mesh_.number_of_cells() * data_.number_of_groups() * mesh_.number_of_nodes() * ordinates_.number_of_ordinates(), 0);
            vector<double> phi(mesh_.number_of_cells() * data_.number_of_groups() * mesh_.number_of_nodes(), 1); // average scalar flux
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

                slab_sweep(psi,
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
                    iterations_ = it + 1;
                    phi_ = phi;
                }
            }
            // calculate_leakage(psi,
            //                   leakage);
        }
        else if (mesh_.geometry("spherical"))
        {
            // temporary variables
            vector<double> psi(mesh_.number_of_cells() * data_.number_of_groups() * mesh_.number_of_nodes() * ordinates_.number_of_ordinates(), 0);
            vector<double> phi(mesh_.number_of_cells() * data_.number_of_groups() * mesh_.number_of_nodes(), 1); // average scalar flux
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
                
                spherical_sweep(psi,
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
                    iterations_ = it + 1;
                    phi_ = phi;
                }
            }

            // print_angular_flux(psi);
            // calculate_leakage(psi,
            //                   leakage);
        }
    }

    void Sn_Transport::
    calculate_leakage(vector<double> &psi,
                      vector<double> &leakage)
    {
        unsigned il = 0;
        unsigned ir = mesh_.number_of_cells() - 1;

        leakage.assign(leakage.size(), 0.0);
        
        for (unsigned on = 0; on < ordinates_.number_of_ordinates() / 2; ++on)
        {
            unsigned op = ordinates_.number_of_ordinates() - on - 1;
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                unsigned kl = 0 + mesh_.number_of_nodes() * (g + data_.number_of_groups() * (op + ordinates_.number_of_ordinates() * il));
                unsigned kr = 1 + mesh_.number_of_nodes() * (g + data_.number_of_groups() * (on + ordinates_.number_of_ordinates() * ir));
                
                leakage[0] += ordinates_.ordinates(op) * ordinates_.weights(op) * psi[kl];
                leakage[1] += ordinates_.ordinates(on) * ordinates_.weights(on) * psi[kr];
            }
        }
    }
}
