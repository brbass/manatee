#include "FEM_Transport.hh"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseSolver.h>

#include "Data.hh"
#include "Mesh.hh"

namespace transport_ns
{
    using namespace std;
    using namespace data_ns;
    using namespace mesh_ns;
    
    FEM_Transport::
    FEM_Transport(Data &data,
                  Mesh &mesh,
                  unsigned max_iterations,
                  double tolerance)
    :
        data_(data),
        mesh_(mesh),
        max_iterations_(max_iterations),
        tolerance_(tolerance)
    {
        data_.check();
        mesh_.check();
    }
    
    void FEM_Transport::
    calculate_source()
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
                        sum += (data_.chi(i, g) * data_.nu_sigma_f(i, g2) + data_.sigma_s(i, g2, g)) * phi(i, g, n);
                    }
                    
                    q(i, g, n) = (data_.internal_source(i, g) + sum) / 2;
                }
            }
        }
    }
    
    void FEM_Transport::
    spherical_sweep()
    {
        // boundary condition, r=R
        if (data_.boundary_condition(1).compare("reflected") == 0)
        {
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                unsigned m = 0;
                unsigned p = 1;
                unsigned i = mesh_.number_of_cells() - 1;
                unsigned n = 1;
                
                for (unsigned o = 0; o < number_of_angular_nodes_; ++o)
                {
                    psi_boundary_sources(g, o, m) = psi(i, g, n, number_of_angular_nodes_ - o - 1, p);
                }
            }
        }
        else if (data_.boundary_condition(1).compare("vacuum") == 0)
        {
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                unsigned m = 0;
                
                for (unsigned o = 0; o < number_of_angular_nodes_; ++o)
                {
                    psi_boundary_sources(g, o, m) = 0;
                }
            }
        }
        else
        {
            cerr << "boundary condition \"" << data_.boundary_condition(1) << "\" not available" << endl;
        }

        // sweep inward
        
        sweep_inward();
        
        // boundary condition, r=0
        if (data_.boundary_condition(0).compare("reflected") == 0)
        {
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                unsigned m = 0;
                unsigned p = 1;
                unsigned i = 0;
                unsigned n = 0;
                
                for (unsigned o = 0; o < number_of_angular_nodes_; ++o)
                {
                    psi_boundary_sources(g, o, p) = psi(i, g, n, number_of_angular_nodes_ - o - 1, m);
                }
            }
        }
        else
        {
            cerr << "boundary condition \"" << data_.boundary_condition(0) << "\" not available" << endl;
        }
        
        sweep_outward();

        // print_angular_flux();
    }

    void FEM_Transport::
    sweep_inward()
    {
        matrix_.assign(total_nodes_ * total_nodes_, 0);
        lhs_.assign(total_nodes_, 0);
        rhs_.assign(total_nodes_, 0);
        
        unsigned m = 0;
        unsigned p = 1;
        
        // first cell
        for (unsigned g = 0; g < data_.number_of_groups(); ++g)
        {
            unsigned i = mesh_.number_of_cells() - 1;
            
            for (unsigned nb = 0; nb < mesh_.number_of_nodes(); ++nb)
            {
                for (unsigned ob = 0; ob < number_of_angular_nodes_; ++ob)
                {
                    unsigned be = ob + number_of_angular_nodes_ * nb;
                        
                    double sum = 0;
                        
                    for (unsigned na = 0; na < mesh_.number_of_nodes(); ++na)
                    {
                        for (unsigned oa = 0; oa < number_of_angular_nodes_; ++oa)
                        {
                            unsigned al = oa + number_of_angular_nodes_ * na;
                            
                            matrix(al, be) = get_k_l(i, al, be, m) + get_l(i, al, be, m) + get_m(i, al, be) * data_.sigma_t(i, g);
                            
                            sum += get_m(i, al, be) * q(i, g, na) / 2;
                            sum -= get_j(i, al, be, m) * psi(i, g, na, oa, p);
                            sum -= get_k_r(i, al, be, m) * psi_boundary_sources(g, oa, m);
                        }
                    }
                    
                    rhs(be) = sum;
                }
            }

            epetra_solve(matrix_, lhs_, rhs_, total_nodes_);
                
            for (unsigned n = 0; n < mesh_.number_of_nodes(); ++n)
            {
                for (unsigned o = 0; o < number_of_angular_nodes_; ++o)
                {
                    psi(i, g, n, o, m) = lhs(n, o);
                }
            }
        }                

        // remaining cells
        for (int i = mesh_.number_of_cells() - 2; i >=0; --i)
        {
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                for (unsigned nb = 0; nb < mesh_.number_of_nodes(); ++nb)
                {
                    for (unsigned ob = 0; ob < number_of_angular_nodes_; ++ob)
                    {
                        unsigned be = ob + number_of_angular_nodes_ * nb;
                        
                        double sum = 0;
                        
                        for (unsigned na = 0; na < mesh_.number_of_nodes(); ++na)
                        {
                            for (unsigned oa = 0; oa < number_of_angular_nodes_; ++oa)
                            {
                                unsigned al = oa + number_of_angular_nodes_ * na;

                                matrix(al, be) = get_k_l(i, al, be, m) + get_l(i, al, be, m) + get_m(i, al, be) * data_.sigma_t(i, g);
                                
                                sum += get_m(i, al, be) * q(i, g, na) / 2;
                                sum -= get_j(i, al, be, m) * psi(i, g, na, oa, p);
                                sum -= get_k_r(i, al, be, m) * psi(i+1, g, 0, oa, m);
                            }
                        }
                        
                        rhs(be) = sum;
                    }
                }

                cout << "*************************************************" << endl;
                cout << "i=" << i << " ri=" << mesh_.cell_center_position(i) << " dri=" << mesh_.cell_length(i) << " s_t=" << data_.sigma_t(i, g) << endl << endl;
                
                epetra_solve(matrix_, lhs_, rhs_, total_nodes_);
                
                for (unsigned n = 0; n < mesh_.number_of_nodes(); ++n)
                {
                    for (unsigned o = 0; o < number_of_angular_nodes_; ++o)
                    {
                        psi(i, g, n, o, m) = lhs(n, o);
                    }
                }
            }
        }
        
    }

    void FEM_Transport::
    sweep_outward()
    {
        matrix_.assign(total_nodes_ * total_nodes_, 0);
        lhs_.assign(total_nodes_, 0);
        rhs_.assign(total_nodes_, 0);

        unsigned m = 0;
        unsigned p = 1;

        // first cell
        for (unsigned g = 0; g < data_.number_of_groups(); ++g)
        {
            unsigned i = 0;
            for (unsigned nb = 0; nb < mesh_.number_of_nodes(); ++nb)
            {
                for (unsigned ob = 0; ob < number_of_angular_nodes_; ++ob)
                {
                    unsigned be = ob + number_of_angular_nodes_ * nb;
                        
                    double sum = 0;
                        
                    for (unsigned oa = 0; oa < number_of_angular_nodes_; ++oa)
                    {
                        for (unsigned na = 0; na < mesh_.number_of_nodes(); ++na)
                        {
                            unsigned al = oa + number_of_angular_nodes_ * na;

                            matrix(al, be) = get_k_r(i, al, be, p) + get_l(i, al, be, p) + get_m(i, al, be) * data_.sigma_t(i, g);
                                
                            sum += get_m(i, al, be) * q(i, g, na) / 2;
                            sum -= get_j(i, al, be, p) * psi(i, g, na, oa, m);
                            sum -= get_k_l(i, al, be, p) * psi_boundary_sources(g, oa, p);
                        }
                    }
                    
                    rhs(be) = sum;
                }
            }

            epetra_solve(matrix_, lhs_, rhs_, total_nodes_);
                
            for (unsigned n = 0; n < mesh_.number_of_nodes(); ++n)
            {
                for (unsigned o = 0; o < number_of_angular_nodes_; ++o)
                {
                    psi(i, g, n, o, p) = lhs(n, o);
                }
            }
        }                

        // remaining cells
        for (unsigned i = 1; i < mesh_.number_of_cells(); ++i)
        {
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                for (unsigned nb = 0; nb < mesh_.number_of_nodes(); ++nb)
                {
                    for (unsigned ob = 0; ob < number_of_angular_nodes_; ++ob)
                    {
                        unsigned be = ob + number_of_angular_nodes_ * nb;
                        
                        double sum = 0;
                        
                        for (unsigned na = 0; na < mesh_.number_of_nodes(); ++na)
                        {
                            for (unsigned oa = 0; oa < number_of_angular_nodes_; ++oa)
                            {
                                unsigned al = oa + number_of_angular_nodes_ * na;

                                matrix(al, be) = get_k_r(i, al, be, p) + get_l(i, al, be, p) + get_m(i, al, be) * data_.sigma_t(i, g);
                                
                                sum += get_m(i, al, be) * q(i, g, na) / 2;
                                sum -= get_j(i, al, be, p) * psi(i, g, na, oa, m);
                                sum -= get_k_l(i, al, be, p) * psi(i-1, g, mesh_.number_of_nodes() - 1, oa, p);
                            }
                        }
                        
                        rhs(be) = sum;
                    }
                }

                cout << "*************************************************" << endl;
                cout << "i=" << i << " ri=" << mesh_.cell_center_position(i) << " dri=" << mesh_.cell_length(i) << " s_t=" << data_.sigma_t(i, g) << endl << endl;
                
                epetra_solve(matrix_, lhs_, rhs_, total_nodes_);
                
                for (unsigned n = 0; n < mesh_.number_of_nodes(); ++n)
                {
                    for (unsigned o = 0; o < number_of_angular_nodes_; ++o)
                    {
                        psi(i, g, n, o, p) = lhs(n, o);
                    }
                }
            }
        }
    }
    
    void FEM_Transport::
    psi_to_phi()
    {
        for (unsigned i = 0; i < mesh_.number_of_cells(); ++i)
        {
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                for (unsigned n = 0; n < mesh_.number_of_nodes(); ++n)
                {
                    double sum = 0;
                    
                    for (unsigned o = 0; o < number_of_angular_nodes_; ++o)
                    {
                        for (unsigned mp = 0; mp < 2; ++mp)
                        {
                            sum += psi(i, g, n, o, mp);
                        }
                    }
                    
                    phi(i, g, n) = sum;
                }
            }
        }
    }
    
    void FEM_Transport::
    check_convergence()
    {
        error_phi_old_ = error_phi_;
        
        converged_ = true;
        
        for (unsigned i = 0; i < mesh_.number_of_cells(); ++i)
        {
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                for (unsigned n = 0; n < mesh_.number_of_nodes(); ++n)
                {
                    error_phi(i, g, n) = abs((phi(i, g, n) - phi_old(i, g, n)) / phi(i, g, n));
                    
                    if (error_phi(i, g, n) > tolerance_ * (1 - error_phi(i, g, n) / error_phi_old(i, g, n)))
                    {
                        converged_ = false;
                    }
                }
            }
        }
    }
    
    void FEM_Transport::
    analytic_solve(vector<double> &matrix, vector<double> &lhs, vector<double> &rhs)
    {
        lhs[0] = (matrix[3] * rhs[0] - matrix[1] * rhs[1]) / (matrix[0] * matrix[3] - matrix[1] * matrix[2]);
        lhs[1] = (matrix[2] * rhs[0] - matrix[0] * rhs[1]) / (matrix[1] * matrix[2] - matrix[0] * matrix[3]);
    }
    
    void FEM_Transport::
    epetra_solve(vector<double> &matrix, vector<double> &lhs, vector<double> &rhs, int size)
    {
        Epetra_SerialDenseMatrix epetra_matrix (View, &matrix[0], size, size, size);
        Epetra_SerialDenseVector epetra_lhs (size);
        Epetra_SerialDenseVector epetra_rhs (View, &rhs[0], size);
        
        Epetra_SerialDenseSolver solver;
        
        cout << epetra_matrix << endl;
        cout << epetra_rhs << endl;
        
        solver.FactorWithEquilibration(true);
        solver.SetMatrix(epetra_matrix);
        solver.SetVectors(epetra_lhs, epetra_rhs);
        solver.Factor();
        solver.Solve();
        
        lhs.assign(solver.X(), solver.X() + size);
        
        cout << "Solution" << endl;
        for (unsigned i = 0; i < 4; ++i)
        {
            cout << lhs[i] << "    ";
        }
        cout << endl << endl;
    }
    
    void FEM_Transport::
    solve()
    {
        using namespace std;
        
        psi_.resize(mesh_.number_of_cells() * data_.number_of_groups() * mesh_.number_of_nodes() * number_of_angular_nodes_ * 2, 0);
        phi_.resize(mesh_.number_of_cells() * data_.number_of_groups() * mesh_.number_of_nodes(), 1);
        phi_old_.resize(mesh_.number_of_cells() * data_.number_of_groups() * mesh_.number_of_nodes(), 0);
        q_.resize(mesh_.number_of_cells() * data_.number_of_groups() * mesh_.number_of_nodes(), 0); // total source, including fission and scattering
        error_phi_.resize(mesh_.number_of_cells() * data_.number_of_groups() * mesh_.number_of_nodes(), 1);
        error_phi_old_.resize(mesh_.number_of_cells() * data_.number_of_groups() * mesh_.number_of_nodes(), 1);
        psi_boundary_sources_.resize(2 * data_.number_of_groups() * number_of_angular_nodes_, 0);

        // for (unsigned mp = 0; mp < 2; ++mp)
        // {
        //     for (int i = mesh_.number_of_cells() - 2; i >=0; --i)
        //     {
        //         for (unsigned g = 0; g < data_.number_of_groups(); ++g)
        //         {
        //             cout << "*************************************************" << endl;
        //             cout << "mp=" << mp << " i=" << i << " ri=" << mesh_.cell_center_position(i) << " dri=" << mesh_.cell_length(i) << " s_t=" << data_.sigma_t(i, g) << endl << endl;
                
        //             for (unsigned nb = 0; nb < mesh_.number_of_nodes(); ++nb)
        //             {
        //                 for (unsigned ob = 0; ob < number_of_angular_nodes_; ++ob)
        //                 {
        //                     unsigned be = ob + number_of_angular_nodes_ * nb;
                        
        //                     double sum = 0;
                        
        //                     for (unsigned na = 0; na < mesh_.number_of_nodes(); ++na)
        //                     {
        //                         for (unsigned oa = 0; oa < number_of_angular_nodes_; ++oa)
        //                         {
        //                             unsigned al = oa + number_of_angular_nodes_ * na;
                                
        //                             cout << "al=" << al << " be=" << be << " j=" << get_j(i, al, be, mp) << " kr=" << get_k_r(i, al, be, mp) << " kl=" << get_k_l(i, al, be, mp) << " l=" << get_l(i, al, be, mp) << " m=" <<  get_m(i, al, be) << endl;
        //                         }
        //                     }
        //                 }
        //             }

        //         }
        //     }
        // }
        
        // begin iterations
        for (unsigned it = 0; it < max_iterations_; ++it)
        {
            calculate_source();
            
            spherical_sweep();
                
            phi_old_ = phi_;
                
            psi_to_phi();
            
            check_convergence();
            
            if (converged_)
            {
                iterations_ = it + 1;
                
                break;
            }
            else if (it == max_iterations_ - 1)
            {
                iterations_ = it + 1;
            }
        }
        
        // print_scalar_flux();
        // print_angular_flux();
        // calculate_leakage();
    }

    void FEM_Transport::
    solve_eigenvalue()
    {
        using namespace std;
        
        psi_.resize(mesh_.number_of_cells() * data_.number_of_groups() * mesh_.number_of_nodes() * number_of_angular_nodes_ * 2, 0);
        phi_.resize(mesh_.number_of_cells() * data_.number_of_groups() * mesh_.number_of_nodes(), 1);
        phi_old_.resize(mesh_.number_of_cells() * data_.number_of_groups() * mesh_.number_of_nodes(), 0);
        q_.resize(mesh_.number_of_cells() * data_.number_of_groups() * mesh_.number_of_nodes(), 0); // total source, including fission and scattering
        error_phi_.resize(mesh_.number_of_cells() * data_.number_of_groups() * mesh_.number_of_nodes(), 1);
        error_phi_old_.resize(mesh_.number_of_cells() * data_.number_of_groups() * mesh_.number_of_nodes(), 1);
        psi_boundary_sources_.resize(2 * data_.number_of_groups() * number_of_angular_nodes_, 0);
        k_eigenvalue_ = 1.0;
        
        // begin iterations
        for (unsigned it = 0; it < max_iterations_; ++it)
        {
            calculate_source();
            
            spherical_sweep();
                
            phi_old_ = phi_;
            psi_to_phi();

            //k_eigenvalue_old_ = k_eigenvalue_;
            calculate_k();

            normalize_phi();
            
            check_convergence();
            
            if (converged_)
            {
                iterations_ = it + 1;
                
                break;
            }
            else if (it == max_iterations_ - 1)
            {
                iterations_ = it + 1;
            }
        }
        
        // print_scalar_flux();
        // print_angular_flux();
        // calculate_leakage();
    }

    double FEM_Transport::
    get_j(unsigned i, unsigned al, unsigned be, unsigned mp)
    {
        double t1 = - mesh_.cell_length(i) / 96;
        double t2 = 2 * mp_coord[mp] -2 * eta_coord[be] + 2 * eta_coord[al] * (-1 + eta_coord[be] * mp_coord[mp]);
        double t3 = mesh_.cell_length(i) * (eta_coord[al] + eta_coord[be]) + 2 * mesh_.cell_center_position(i) * (3 + xi_coord[al] * xi_coord[be]);

        return t1 * t2 * t3;
           
    }
    
    double FEM_Transport::
    get_k_r(unsigned i, unsigned al, unsigned be, unsigned mp)
    {
        double t1 = 1.0 / 192.0;
        double t2 = eta_coord[al] + eta_coord[be] + 3 * mp_coord[mp] + eta_coord[al] * eta_coord[be] * mp_coord[mp];
        double t3 = pow(mesh_.cell_length(i) + 2 * mesh_.cell_center_position(i), 2);
        double t4 = (xi_coord[al] + 1) * (xi_coord[be] + 1);
        
        return t1 * t2 * t3 * t4;
    }
    
    double FEM_Transport::
    get_k_l(unsigned i, unsigned al, unsigned be, unsigned mp)
    {
        double t1 = - 1.0 / 192.0;
        double t2 = eta_coord[al] + eta_coord[be] + 3 * mp_coord[mp] + eta_coord[al] * eta_coord[be] * mp_coord[mp];
        double t3 = pow(mesh_.cell_length(i) - 2 * mesh_.cell_center_position(i), 2);
        double t4 = (xi_coord[al] - 1) * (xi_coord[be] - 1);
        
        return t1 * t2 * t3 * t4;
    }
    
    double FEM_Transport::
    get_l(unsigned i, unsigned al, unsigned be, unsigned mp)
    {
        double t1 = - 1.0 / 288.0;
        double t2 = eta_coord[al] + eta_coord[be] + 3 * mp_coord[mp] + eta_coord[al] * eta_coord[be] * mp_coord[mp];
        double t3 = (pow(mesh_.cell_length(i), 2) + 12 * pow(mesh_.cell_center_position(i), 2) + 4 * mesh_.cell_length(i) * mesh_.cell_center_position(i) * xi_coord[al]) * xi_coord[be];
        double t4 = mesh_.cell_length(i) / 288.0;
        double t5 = eta_coord[be] * (-8 + 2 * eta_coord[al] * mp_coord[mp]);
        double t6 = mesh_.cell_length(i) * (xi_coord[al] + xi_coord[be]) + 2 * mesh_.cell_center_position(i) * (3 + xi_coord[al] * xi_coord[be]);
        
        return t1 * t2 * t3 + t4 * t5 * t6;
    }
    
    double FEM_Transport::
    get_m(unsigned i, unsigned al, unsigned be)
    {
        double t1 = mesh_.cell_length(i) / 1440.0;
        double t2 = 3 + eta_coord[al] * eta_coord[be];
        double t3 = 20 * mesh_.cell_length(i) * mesh_.cell_center_position(i) * (xi_coord[al] + xi_coord[be]) + 20 * pow(mesh_.cell_center_position(i), 2) * (3 + xi_coord[al] * xi_coord[be]) + pow(mesh_.cell_length(i), 2) * (5 + 3 * xi_coord[al] * xi_coord[be]);
        
        return t1 * t2 * t3;
    }

    void FEM_Transport::
    calculate_k()
    {
        double num = 0;
        double den = 0;
        
        for (unsigned i = 0; i < mesh_.number_of_cells(); ++i)
        {
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                for (unsigned n = 0; n < mesh_.number_of_nodes(); ++n)
                {
                    unsigned k1 = n + mesh_.number_of_nodes() * (g + data_.number_of_groups() * i);
                    
                    num += mesh_.cell_length(i) * data_.nu_sigma_f(i, g) * phi_[k1];
                    den += mesh_.cell_length(i) * data_.nu_sigma_f(i, g) * phi_old_[k1];
                }
            }
        }
        
        k_eigenvalue_ = /*k_eigenvalue_old_ **/ num / den;
    }

    void FEM_Transport::
    normalize_phi()
    {
        double sum = 0;
        
        for (unsigned i = 0; i < mesh_.number_of_cells(); ++i)
        {
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                for (unsigned n = 0; n < mesh_.number_of_nodes(); ++n)
                {
                    sum += mesh_.cell_length(i) * phi(i, g, n);
                }
            }
        }

        for (unsigned i = 0; i < mesh_.number_of_cells(); ++i)
        {
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                for (unsigned n = 0; n < mesh_.number_of_nodes(); ++n)
                {
                    phi(i, g, n) = abs(phi(i, g, n)) / sum;
                }
            }
        }
    }

}
