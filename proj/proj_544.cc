#include <iostream>
#include <cmath>
#include <string>

#include <Teuchos_Time.hpp>

#include <Monte_Carlo.hh>
#include <SP1_Transport.hh>

int test_weighted_mc()
{
    using namespace std;
    using namespace monte_carlo_ns;
    using namespace transport_ns;
    
    // initialize mutual variables
    unsigned number_of_cells = 30;
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
    
    SP1_Transport sp1_transport(number_of_cells,
                                number_of_groups,
                                number_of_scattering_moments,
                                side_length,
                                response, // use cross section in detector as response
                                boundary_sources,
                                sigma_t,
                                sigma_s,
                                nu_sigma_f,
                                chi,
                                boundary_conditions,
                                problem_type);
    
    sp1_transport.solve();
    
    sp1_transport.print_scalar_flux();
    
    vector<double> phi_adjoint = sp1_transport.get_scalar_flux();
    
    // solve using monte carlo and weight windows
    boundary_conditions.assign(2, "reflected");
    unsigned number_of_histories = 1e6;
    bool implicit_capture = true;
    
    vector<double> nu(number_of_cells * number_of_groups, 0.0*1.70);
    vector<double> sigma_f(number_of_cells * number_of_groups, 0.0*0.054628);
    
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
    
    monte_carlo.initialize_weight_windows(phi_adjoint);
    
    monte_carlo.solve();
    
    monte_carlo.print_scalar_flux();
    
    return 0;
}


int main(int argc, char *argv[])
{
    test_weighted_mc();
}
