#include <iomanip>
#include <iostream>
#include <string>

#include <Teuchos_Time.hpp>

#include "Sn_Transport.hh"

/*
Critical bare homogeneous U-D2O system with isotropic scattering
    v = 1.70
    sf = 0.054628
    ss = 0.464338
    st = 0.54628
    R = 22.017156
*/
int u_d2o_sys()
{
    using namespace std;
    using namespace transport_ns;

    vector<unsigned> numbers_of_cells = {};
    for (unsigned i = 1; i < 101; ++i)
    {
        numbers_of_cells.push_back(i);
    }
    
    vector<unsigned> numbers_of_ordinates = {};
    for (unsigned o = 1; o < 33; ++o)
    {
        numbers_of_ordinates.push_back(2 * o);
    }
    
    vector<string> boundary_conditions(2, "reflected");
    boundary_conditions[1] = "vacuum";
    string geometry = "spherical";

    cout << endl << "U-D2O System" << endl;
    
    for (unsigned i = 0; i < numbers_of_cells.size(); ++i)
    {
        for (unsigned o = 0; o < numbers_of_ordinates.size(); ++o)
        {
            unsigned number_of_cells = numbers_of_cells[i];
            unsigned number_of_groups = 1;
            unsigned number_of_scattering_moments = 1;
            unsigned number_of_ordinates = numbers_of_ordinates[o];
    
            unsigned max_iterations = 10000;
            double tolerance = 1e-10;
    
            double side_length = 22.017156;
            vector<double> internal_source(number_of_cells * number_of_groups, 0.0);
            vector<double> boundary_sources(number_of_groups * number_of_ordinates * 2, 0.0);
            vector<double> sigma_t(number_of_cells * number_of_groups, 0.54628);
            vector<double> sigma_s(number_of_cells * number_of_groups * number_of_groups * number_of_scattering_moments, 0.464338);
            vector<double> nu_sigma_f(number_of_cells * number_of_groups, 1.70*0.054628);
            vector<double> chi(number_of_cells * number_of_groups, 1.0);
    
            Sn_Transport sn_transport(number_of_cells,
                                      number_of_groups,
                                      number_of_scattering_moments,
                                      number_of_ordinates,
                                      side_length,
                                      internal_source,
                                      boundary_sources,
                                      sigma_t,
                                      sigma_s,
                                      nu_sigma_f,
                                      chi,
                                      boundary_conditions,
                                      geometry);

            Teuchos::Time timer("proj");
            timer.start();
            sn_transport.solve_eigenvalue();
            double t = timer.stop();
            
            double k = sn_transport.get_eigenvalue();
            double it = sn_transport.get_iterations();
            
            cout << number_of_cells << "\t" << number_of_ordinates << "\t" << it << "\t" << setw(10) << setprecision(8) << t << "\t" << setprecision(25) << k << endl;
        }
    }

    return 0;
}

// Test two independent, identical groups.
int u_d2o_test()
{
    using namespace std;
    using namespace transport_ns;

    vector<unsigned> numbers_of_cells = {10};
    // for (unsigned i = 1; i < 21; ++i)
    // {
    //     numbers_of_cells.push_back(i);
    // }
    
    vector<unsigned> numbers_of_ordinates = {16};
    // for (unsigned o = 1; o < 9; ++o)
    // {
    //     numbers_of_ordinates.push_back(2 * o);
    // }

    vector<string> boundary_conditions(2, "reflected");
    boundary_conditions[1] = "vacuum";
    string geometry = "spherical";

    cout << endl << "U-D2O Two-Group Test" << endl;
    
    for (unsigned i = 0; i < numbers_of_cells.size(); ++i)
    {
        for (unsigned o = 0; o < numbers_of_ordinates.size(); ++o)
        {
            unsigned number_of_cells = numbers_of_cells[i];
            unsigned number_of_groups = 2;
            unsigned number_of_scattering_moments = 1;
            unsigned number_of_ordinates = numbers_of_ordinates[o];
    
            unsigned max_iterations = 10000;
            double tolerance = 1e-10;
    
            double side_length = 22.017156;
            vector<double> internal_source(number_of_cells * number_of_groups, 0.0);
            vector<double> boundary_sources(number_of_groups * number_of_ordinates * 2, 0.0);
            vector<double> sigma_t(number_of_cells * number_of_groups, 0.54628);
            vector<double> sigma_s(number_of_cells * number_of_groups * number_of_groups * number_of_scattering_moments, 0.464338 / number_of_groups);
            vector<double> nu_sigma_f(number_of_cells * number_of_groups, 1.70*0.054628);
            vector<double> chi(number_of_cells * number_of_groups, 1.0 / number_of_groups);
    
            Sn_Transport sn_transport(number_of_cells,
                                      number_of_groups,
                                      number_of_scattering_moments,
                                      number_of_ordinates,
                                      side_length,
                                      internal_source,
                                      boundary_sources,
                                      sigma_t,
                                      sigma_s,
                                      nu_sigma_f,
                                      chi,
                                      boundary_conditions,
                                      geometry);
    
            Teuchos::Time timer("proj");
            timer.start();
            sn_transport.solve_eigenvalue();
            double t = timer.stop();
            
            double k = sn_transport.get_eigenvalue();
            double it = sn_transport.get_iterations();
            
            cout << number_of_cells << "\t" << number_of_ordinates << "\t" << it << "\t" << setw(10) << setprecision(8) << t << "\t" << setprecision(25) << k << endl;
        }
    }

    return 0;
}

/*
Two-group model of a critical bare, homogeneous system made of 93% enriched uranium
    v = 2.5, 2.5
    sf = 0.0010484, 0.050632
    ss1g = 0.62568, 0.029227
    ss2g = 0.0, 2.44383
    st = 0.65696, 2.52025
    chi = 1.0, 0.0
    R = 16.049836
*/
int uranium_sys()
{
    using namespace std;
    using namespace transport_ns;
    
    vector<unsigned> numbers_of_cells = {};
    for (unsigned i = 1; i < 101; ++i)
    {
        numbers_of_cells.push_back(i);
    }
    
    vector<unsigned> numbers_of_ordinates = {};
    for (unsigned o = 1; o < 33; ++o)
    {
        numbers_of_ordinates.push_back(2 * o);
    }
    
    unsigned number_of_groups = 2;
    unsigned number_of_scattering_moments = 1;
    unsigned max_iterations = 10000;
    double tolerance = 1e-10;
    double side_length = 16.049836;
    
    vector<string> boundary_conditions(2, "reflected");
    boundary_conditions[1] = "vacuum";
    string geometry = "spherical";
    
    cout << endl << "Uranium System" << endl;
    
    for (unsigned i1 = 0; i1 < numbers_of_cells.size(); ++i1)
    {
        for (unsigned o1 = 0; o1 < numbers_of_ordinates.size(); ++o1)
        {
            unsigned number_of_cells = numbers_of_cells[i1];
            unsigned number_of_ordinates = numbers_of_ordinates[o1];
            
            vector<double> internal_source(number_of_cells * number_of_groups, 0.0);
            vector<double> boundary_sources(number_of_groups * number_of_ordinates * 2, 0.0);
            vector<double> sigma_t(number_of_cells * number_of_groups, 0);
            vector<double> sigma_s(number_of_cells * number_of_groups * number_of_groups * number_of_scattering_moments, 0);
            vector<double> nu_sigma_f(number_of_cells * number_of_groups, 0);
            vector<double> chi(number_of_cells * number_of_groups, 0.0);
            for (unsigned i = 0; i < number_of_cells; ++i)
            {
                unsigned k0 = 0 + number_of_groups * i;
                unsigned k1 = 1 + number_of_groups * i;
                unsigned k00 = 0 + number_of_groups * (0 + number_of_groups * i);
                unsigned k01 = 0 + number_of_groups * (1 + number_of_groups * i);
                unsigned k10 = 1 + number_of_groups * (0 + number_of_groups * i);
                unsigned k11 = 1 + number_of_groups * (1 + number_of_groups * i);
                
                sigma_t[k0] = 0.65696;
                sigma_t[k1] = 2.52025;
                
                nu_sigma_f[k0] = 2.5 * 0.0010484;
                nu_sigma_f[k1] = 2.5 * 0.050632;
                
                internal_source[k0] = 0.0;
                internal_source[k1] = 0.0;
                
                chi[k0] = 1.0;
                chi[k1] = 0.0;
                
                sigma_s[k00] = 0.62568;
                sigma_s[k01] = 0.029227;
                sigma_s[k10] = 0.0;
                sigma_s[k11] = 2.44383;
            }
            
            Sn_Transport sn_transport(number_of_cells,
                                      number_of_groups,
                                      number_of_scattering_moments,
                                      number_of_ordinates,
                                      side_length,
                                      internal_source,
                                      boundary_sources,
                                      sigma_t,
                                      sigma_s,
                                      nu_sigma_f,
                                      chi,
                                      boundary_conditions,
                                      geometry);
            
            Teuchos::Time timer("proj");
            timer.start();
            sn_transport.solve_eigenvalue();
            double t = timer.stop();
            
            double k = sn_transport.get_eigenvalue();
            double it = sn_transport.get_iterations();
            
            cout << number_of_cells << "\t" << number_of_ordinates << "\t" << it << "\t" << setw(10) << setprecision(8) << t << "\t" << setprecision(25) << k << endl;
        }
    }
    
    return 0;
}

/*
One-group model of a critical, water reflected U-235 fuelled system
  u-235:
    v = 2.679198
    sf = 0.065280
    ss = 0.248064
    st = 0.36240
    R = 0, 6.12745
  h2o:
    v = 0.0
    sf = 0.0
    ss = 0.293760
    st = 0.32640
    R = 6.12745, 15.318626
*/
int u_235_sys()
{
    using namespace std;
    using namespace transport_ns;

    vector<unsigned> numbers_of_cells = {};
    for (unsigned i = 1; i < 21; ++i)
    {
        numbers_of_cells.push_back(i * 5);
    }
    
    vector<unsigned> numbers_of_ordinates = {};
    for (unsigned o = 1; o < 33; ++o)
    {
        numbers_of_ordinates.push_back(2 * o);
    }

    vector<string> boundary_conditions(2, "reflected");
    boundary_conditions[1] = "vacuum";
    string geometry = "spherical";
    
    cout << endl << "U-235 System" << endl;
    
    for (unsigned i1 = 0; i1 < numbers_of_cells.size(); ++i1)
    {
        for (unsigned o1 = 0; o1 < numbers_of_ordinates.size(); ++o1)
        {
            unsigned number_of_cells = numbers_of_cells[i1];

            if (number_of_cells % 5 != 0)
            {
                cerr << "only powers of five supported" << endl;
            }

            unsigned mid_cell = number_of_cells / 5 * 2;

            unsigned number_of_groups = 1;
            unsigned number_of_scattering_moments = 1;
            unsigned number_of_ordinates = numbers_of_ordinates[o1];
            
            unsigned max_iterations = 10000;
            double tolerance = 1e-10;
            
            double side_length = 15.318626;
            vector<double> internal_source(number_of_cells * number_of_groups, 0.0);
            vector<double> boundary_sources(number_of_groups * number_of_ordinates * 2, 0.0);
            vector<double> sigma_t(number_of_cells * number_of_groups, 0.0);
            vector<double> sigma_s(number_of_cells * number_of_groups * number_of_groups * number_of_scattering_moments, 0.0);
            vector<double> nu_sigma_f(number_of_cells * number_of_groups, 0.0);
            vector<double> chi(number_of_cells * number_of_groups, 0.0);

            for (unsigned i = 0; i < number_of_cells; ++i)
            {
                unsigned k0 = 0 + number_of_groups * i;
                unsigned k00 = 0 + number_of_groups * (0 + number_of_groups * i);

                if (i < mid_cell)
                {
                    sigma_t[k0] = 0.32640;
                    nu_sigma_f[k0] = 2.679198 * 0.065280;
                    sigma_s[k00] = 0.248064;
                    chi[k0] = 1;
                }
                else
                {
                    sigma_t[k0] = 0.32640;
                    nu_sigma_f[k0] = 0;
                    sigma_s[k00] = 0.293760;
                    chi[k0] = 0;
                }
            }
            
            Sn_Transport sn_transport(number_of_cells,
                                      number_of_groups,
                                      number_of_scattering_moments,
                                      number_of_ordinates,
                                      side_length,
                                      internal_source,
                                      boundary_sources,
                                      sigma_t,
                                      sigma_s,
                                      nu_sigma_f,
                                      chi,
                                      boundary_conditions,
                                      geometry);
    
            Teuchos::Time timer("proj");
            timer.start();
            sn_transport.solve_eigenvalue();
            double t = timer.stop();
            
            double k = sn_transport.get_eigenvalue();
            double it = sn_transport.get_iterations();
            
            cout << number_of_cells << "\t" << number_of_ordinates << "\t" << it << "\t" << setw(10) << setprecision(8) << t << "\t" << setprecision(25) << k << endl;
        }
    }

    return 0;
}

int main(int argc, char *argv[])
{
    unsigned checksum = 0;

    // u_d2o_sys();
    //u_d2o_test();
    // uranium_sys();
    u_235_sys();
    
    return checksum;
}
