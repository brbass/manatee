#include <iostream>
#include <cmath>
#include <string>

#include <Teuchos_Time.hpp>

#include <Monte_Carlo.hh>
#include <Sn_Transport.hh>
#include <SP1_Transport.hh>

int thick_slab(unsigned number_of_histories)
{
    using namespace std;
    using namespace monte_carlo_ns;
    using namespace transport_ns;
    
    // initialize mutual variables
    unsigned number_of_cells = 300;
    unsigned number_of_groups = 2;
    unsigned number_of_scattering_moments = 2;
    
    double side_length = 10;
    vector<double> internal_source(number_of_cells * number_of_groups, 0.0);
    vector<double> response(number_of_cells * number_of_groups, 0.0);
    vector<double> boundary_sources(number_of_groups * 2, 0.0);
    vector<double> sigma_t(number_of_cells * number_of_groups, 0.0);
    vector<double> sigma_s(number_of_cells * number_of_groups * number_of_groups * number_of_scattering_moments, 0.0);
    vector<double> chi(number_of_cells * number_of_groups, 0.0);
    
    unsigned number_of_regions = 3;
    
    // vector<double> sigma_t0 = {0.1, 10.0, 1.0};
    // vector<double> sigma_t1 = {0.1, 10.0, 1.0};

    // vector<double> sigma_s00 = {0.04, 4.0, 0.4};
    // vector<double> sigma_s01 = {0.04, 4.0, 0.4};
    // vector<double> sigma_s10 = {0.0, 0.0, 0.0};
    // vector<double> sigma_s11 = {0.08, 8.0, 0.0};

    vector<double> sigma_t0 = {0.1, 10.0, 1.0};
    vector<double> sigma_t1 = {0.1, 10.0, 1.0};
    
    vector<double> sigma_s00 = {0.9 * sigma_t0[0], 0.9 * sigma_t0[1], 0.45 * sigma_t0[2]};
    vector<double> sigma_s01 = {0.1 * sigma_t0[0], 0.1 * sigma_t0[1], 0.45 * sigma_t0[2]};
    vector<double> sigma_s10 = {0.0, 0.0, 0.0};
    vector<double> sigma_s11 = {0.9 * sigma_t1[0], 0.9 * sigma_t1[1], 0.0};
    
    vector<double> internal_source0 = {1.0, 0.0, 0.0};
    vector<double> internal_source1 = {0.0, 0.0, 0.0};

    vector<double> response0 = {0.0, 0.0, 0.0};
    vector<double> response1 = {0.0, 0.0, 1.0};
    
    for (unsigned i = 0; i < number_of_cells; ++i)
    {
        unsigned r = floor(1.0 * number_of_regions * i / number_of_cells);
        
        unsigned k0 = 0 + number_of_groups * i;
        unsigned k1 = 1 + number_of_groups * i;
        
        sigma_t[k0] = sigma_t0[r];
        sigma_t[k1] = sigma_t1[r];
        
        internal_source[k0] = internal_source0[r];
        internal_source[k1] = internal_source1[r]; 
        response[k0] = response0[r];
        response[k1] = response1[r];
        
        unsigned k00 = 0 + number_of_groups * (0 + number_of_groups * i);
        unsigned k01 = 0 + number_of_groups * (1 + number_of_groups * i);
        unsigned k10 = 1 + number_of_groups * (0 + number_of_groups * i);
        unsigned k11 = 1 + number_of_groups * (1 + number_of_groups * i);
        unsigned k00m = 0 + number_of_groups * (0 + number_of_groups * (i + number_of_cells * 1));
        unsigned k01m = 0 + number_of_groups * (1 + number_of_groups * (i + number_of_cells * 1));
        unsigned k10m = 1 + number_of_groups * (0 + number_of_groups * (i + number_of_cells * 1));
        unsigned k11m = 1 + number_of_groups * (1 + number_of_groups * (i + number_of_cells * 1));
        
        sigma_s[k00] = sigma_s00[r];
        sigma_s[k01] = sigma_s01[r];
        sigma_s[k10] = sigma_s10[r];
        sigma_s[k11] = sigma_s11[r];
        sigma_s[k00m] = sigma_s[k00] / 4;
        sigma_s[k01m] = sigma_s[k01] / 4;
        sigma_s[k10m] = sigma_s[k10] / 4;
        sigma_s[k11m] = sigma_s[k11] / 4;
    }
    
    // solve sp1 equation for adjoint
    vector<double> nu_sigma_f(number_of_cells * number_of_groups, 0.0);
    vector<string> boundary_conditions(2, "marshak"); // dirichlet, marshak, robin
    string problem_type = "adjoint";
    
    SP1_Transport sp1_transport_ad(number_of_cells,
                                   number_of_groups,
                                   number_of_scattering_moments,
                                   side_length,
                                   response, // use area to which neutrons should be pushed as response
                                   boundary_sources,
                                   sigma_t,
                                   sigma_s,
                                   nu_sigma_f,
                                   chi,
                                   boundary_conditions,
                                   problem_type);
    
    sp1_transport_ad.solve();
    
    sp1_transport_ad.print_scalar_flux();
    
    vector<double> phi_adjoint = sp1_transport_ad.get_scalar_flux();

    // solve sp1 equation for forward flux
    problem_type = "forward";
    
    SP1_Transport sp1_transport_fo(number_of_cells,
                                   number_of_groups,
                                   number_of_scattering_moments,
                                   side_length,
                                   internal_source,
                                   boundary_sources,
                                   sigma_t,
                                   sigma_s,
                                   nu_sigma_f,
                                   chi,
                                   boundary_conditions,
                                   problem_type);
    
    sp1_transport_fo.solve();
    
    sp1_transport_fo.print_scalar_flux();
    
    // solve using monte carlo and weight windows
    boundary_conditions.assign(2, "reflected");
    vector<double> nu(number_of_cells * number_of_groups, 0.0*1.70);
    vector<double> sigma_f(number_of_cells * number_of_groups, 0.0*0.054628);
    
    for (unsigned i = 0; i < 4; ++i)
    {
        bool implicit_capture;
        
        if (i == 1 || i == 3)
        {
            implicit_capture = true;
        }
        else
        {
            implicit_capture = false;
        }
        
        Monte_Carlo monte_carlo(number_of_cells,
                                number_of_groups,
                                number_of_scattering_moments,
                                number_of_histories,
                                side_length,
                                internal_source,
                                boundary_sources,
                                sigma_t,
                                sigma_s,
                                nu,
                                sigma_f,
                                chi,
                                boundary_conditions,
                                implicit_capture);

        if (i == 2 || i == 3)
        {
            monte_carlo.initialize_weight_windows(phi_adjoint);
        }
        
        monte_carlo.solve();
        
        monte_carlo.print_scalar_flux();
        
    }

    return 0;
}

int scattering_slab(unsigned number_of_histories)
{
    using namespace std;
    using namespace monte_carlo_ns;
    using namespace transport_ns;
    
    // initialize mutual variables
    unsigned number_of_cells = 300;
    unsigned number_of_groups = 2;
    unsigned number_of_scattering_moments = 1;
    
    double side_length = 1;
    vector<double> internal_source(number_of_cells * number_of_groups, 0.0);
    vector<double> response(number_of_cells * number_of_groups, 0.0);
    vector<double> boundary_sources(number_of_groups * 2, 0.0);
    vector<double> sigma_t(number_of_cells * number_of_groups, 0.0);
    vector<double> sigma_s(number_of_cells * number_of_groups * number_of_groups * number_of_scattering_moments, 0.0);
    vector<double> nu_sigma_f(number_of_cells * number_of_groups, 0.0);
    vector<double> nu(number_of_cells * number_of_groups, 0.0);
    vector<double> sigma_f(number_of_cells * number_of_groups, 0.0);
    vector<double> chi(number_of_cells * number_of_groups, 0.0);
    string geometry = "slab";
    
    unsigned number_of_regions = 1;
    
    vector<double> sigma_t0 = {10.0};
    vector<double> sigma_t1 = {10.0};
    
    vector<double> sigma_s00 = {0.6 * sigma_t0[0]};
    vector<double> sigma_s01 = {0.39 * sigma_t0[0]};
    vector<double> sigma_s10 = {0.19 * sigma_t1[0]};
    vector<double> sigma_s11 = {0.8};
    
    vector<double> internal_source0 = {1.0};
    vector<double> internal_source1 = {0.0};
    
    vector<double> response0 = {0.0};
    vector<double> response1 = {1.0};
    
    for (unsigned i = 0; i < number_of_cells; ++i)
    {
        unsigned r = floor(1.0 * number_of_regions * i / number_of_cells);
        
        unsigned k0 = 0 + number_of_groups * i;
        unsigned k1 = 1 + number_of_groups * i;
        
        sigma_t[k0] = sigma_t0[r];
        sigma_t[k1] = sigma_t1[r];
        
        internal_source[k0] = internal_source0[r];
        internal_source[k1] = internal_source1[r]; 
        response[k0] = response0[r];
        response[k1] = response1[r];
        
        unsigned k00 = 0 + number_of_groups * (0 + number_of_groups * i);
        unsigned k01 = 0 + number_of_groups * (1 + number_of_groups * i);
        unsigned k10 = 1 + number_of_groups * (0 + number_of_groups * i);
        unsigned k11 = 1 + number_of_groups * (1 + number_of_groups * i);
        unsigned k00m = 0 + number_of_groups * (0 + number_of_groups * (i + number_of_cells * 1));
        unsigned k01m = 0 + number_of_groups * (1 + number_of_groups * (i + number_of_cells * 1));
        unsigned k10m = 1 + number_of_groups * (0 + number_of_groups * (i + number_of_cells * 1));
        unsigned k11m = 1 + number_of_groups * (1 + number_of_groups * (i + number_of_cells * 1));
        
        sigma_s[k00] = sigma_s00[r];
        sigma_s[k01] = sigma_s01[r];
        sigma_s[k10] = sigma_s10[r];
        sigma_s[k11] = sigma_s11[r];
        // sigma_s[k00m] = sigma_s[k00] / 4;
        // sigma_s[k01m] = sigma_s[k01] / 4;
        // sigma_s[k10m] = sigma_s[k10] / 4;
        // sigma_s[k11m] = sigma_s[k11] / 4;
    }
    
    // solve sp1 equation for adjoint
    vector<string> boundary_conditions(2, "robin"); // dirichlet, marshak, robin
    string problem_type = "adjoint";
    
    SP1_Transport sp1_transport_ad(number_of_cells,
                                   number_of_groups,
                                   number_of_scattering_moments,
                                   side_length,
                                   response, // use area to which neutrons should be pushed as response
                                   boundary_sources,
                                   sigma_t,
                                   sigma_s,
                                   nu_sigma_f,
                                   chi,
                                   boundary_conditions,
                                   problem_type);
    
    sp1_transport_ad.solve();
    
    sp1_transport_ad.print_scalar_flux();
    
    vector<double> phi_adjoint = sp1_transport_ad.get_scalar_flux();

    // solve sp1 equation for forward flux
    problem_type = "forward";
    
    SP1_Transport sp1_transport_fo(number_of_cells,
                                   number_of_groups,
                                   number_of_scattering_moments,
                                   side_length,
                                   internal_source,
                                   boundary_sources,
                                   sigma_t,
                                   sigma_s,
                                   nu_sigma_f,
                                   chi,
                                   boundary_conditions,
                                   problem_type);
    
    sp1_transport_fo.solve();
    
    sp1_transport_fo.print_scalar_flux();
    
    // solve using sn transport and then monte carlo with weight windows
    boundary_conditions.assign(2, "vacuum");
    unsigned number_of_ordinates = 256;
    boundary_sources.assign(number_of_groups * number_of_ordinates * 2, 0.0);
    
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

    sn_transport.solve();
    
    sn_transport.print_scalar_flux();
    
    for (unsigned i = 0; i < 4; ++i)
    {
        bool implicit_capture;
        
        if (i == 1 || i == 3)
        {
            implicit_capture = true;
        }
        else
        {
            implicit_capture = false;
        }
        
        Monte_Carlo monte_carlo(number_of_cells,
                                number_of_groups,
                                number_of_scattering_moments,
                                number_of_histories,
                                side_length,
                                internal_source,
                                boundary_sources,
                                sigma_t,
                                sigma_s,
                                nu,
                                sigma_f,
                                chi,
                                boundary_conditions,
                                implicit_capture);

        if (i == 2 || i == 3)
        {
            monte_carlo.initialize_weight_windows(phi_adjoint);
        }
        
        monte_carlo.solve();
        
        monte_carlo.print_scalar_flux();
        
    }

    return 0;
}

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        cerr << "case(0: thick_slab, 1: scattering_slab), number_of_histories" << endl;
        return 1;
    }

    unsigned test_case = atoi(argv[1]);
    unsigned number_of_histories = atoi(argv[2]);

    if (test_case == 0)
    {
        thick_slab(number_of_histories);
    }
    else if (test_case == 1)
    {
        scattering_slab(number_of_histories);
    }
    else
    {
        cerr << "no such test case" << endl;
        return 1;
    }

    return 0;
}
