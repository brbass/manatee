#include <iostream>
#include <string>

#include <Teuchos_Time.hpp>

#include "Data.hh"
#include "Mesh.hh"
#include "Monte_Carlo.hh"
#include "Transport_Model.hh"

int test_homo()
{
    using namespace std;
    using namespace transport_ns;
    
    unsigned number_of_cells = 10;
    unsigned number_of_groups = 1;
    unsigned number_of_scattering_moments = 2;
    
    double side_length = 1;
    vector<unsigned> number_of_cells_per_dimension(1, number_of_cells);
    vector<double> internal_source(number_of_cells * number_of_groups, 1.0);
    vector<double> boundary_sources(number_of_groups * 2, 0.0);
    vector<double> sigma_t(number_of_cells * number_of_groups, 0.5);
    vector<double> sigma_s(number_of_cells * number_of_groups * number_of_groups * number_of_scattering_moments, 0);
    vector<double> nu_sigma_f(number_of_cells * number_of_groups, 0.0);
    vector<double> chi(number_of_cells * number_of_groups, 0.0);
    vector<string> boundary_conditions(2, "marshak"); // dirichlet, marshak, robin
    string problem_type = "forward";
    
    for (unsigned i = 0; i < number_of_cells; ++i)
    {
        for (unsigned gf = 0; gf < number_of_groups; ++gf)
        {
            for (unsigned gt = 0; gt < number_of_groups; ++gt)
            {
                unsigned m = 0;

                unsigned k2 = gf + number_of_groups * (gt + number_of_groups * (i + number_of_cells * m));
                sigma_s[k2] = 0.4 / number_of_groups;
                
                m = 1;

                k2 = gf + number_of_groups * (gt + number_of_groups * (i + number_of_cells * m));
                sigma_s[k2] = 0.2 / number_of_groups;

            }
        }
    }

    Transport_Model transport(number_of_cells,
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
    
    transport.solve();
    
    transport.print_scalar_flux();
    
    return 0;
}

int test_downscatter()
{
    using namespace std;
    using namespace transport_ns;
    
    unsigned number_of_cells = 100;
    unsigned number_of_groups = 2;
    unsigned number_of_scattering_moments = 2;
    unsigned number_of_ordinates = 8;
    
    double side_length = 20;
    vector<unsigned> number_of_cells_per_dimension(1, number_of_cells);
    vector<double> internal_source(number_of_cells * number_of_groups, 0.0);
    vector<double> boundary_sources(number_of_groups * 2, 0.0);
    vector<double> sigma_t(number_of_cells * number_of_groups, 1.0);
    vector<double> sigma_s(number_of_cells * number_of_groups * number_of_groups * number_of_scattering_moments, 0.0);
    vector<double> nu_sigma_f(number_of_cells * number_of_groups, 0.0);
    vector<double> chi(number_of_cells * number_of_groups, 0.0);
    vector<string> boundary_conditions(2, "robin"); // dirichlet, marshak, robin
    
    for (unsigned i = 0; i < number_of_cells; ++i)
    {
        internal_source[0 + number_of_groups * i] = 1;
        
        for (unsigned gf = 0; gf < number_of_groups; ++gf)
        {
            for (unsigned gt = 0; gt < number_of_groups; ++gt)
            {
                unsigned m = 0;
                unsigned k1 = gf + number_of_groups * (gt + number_of_groups * (i + number_of_cells * m));

                m = 1;
                //unsigned k2 = gf + number_of_groups * (gt + number_of_groups * (i + number_of_cells * m));

                if (gf == gt)
                {
                    sigma_s[k1] = 0.4;
                    //sigma_s[k2] = 0.1;
                }
                else if (gf + 1 == gt)
                {
                    sigma_s[k1] = 0.59;
                    //sigma_s[k2] = 0.1;
                }
            }
        }
    }
    
    Transport_Model spn_transport(number_of_cells,
                                  number_of_groups,
                                  number_of_scattering_moments,
                                  side_length,
                                  internal_source,
                                  boundary_sources,
                                  sigma_t,
                                  sigma_s,
                                  nu_sigma_f,
                                  chi,
                                  boundary_conditions);
    
    spn_transport.solve();
    
    spn_transport.print_scalar_flux();

    boundary_conditions.assign(2, "vacuum");
    boundary_sources.assign(number_of_groups * number_of_ordinates * 2, 0.0);
    
    Transport_Model sn_transport(number_of_cells,
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
                             boundary_conditions);

    sn_transport.solve();
    
    sn_transport.print_scalar_flux();
    
    return 0;
}

int test_adjoint()
{
    using namespace std;
    using namespace transport_ns;
    
    unsigned number_of_cells = 5;
    unsigned number_of_groups = 2;
    unsigned number_of_scattering_moments = 2;
    
    double side_length = 1;
    vector<unsigned> number_of_cells_per_dimension(1, number_of_cells);
    vector<double> internal_source(number_of_cells * number_of_groups, 0.0);
    vector<double> boundary_sources(number_of_groups * 2, 0.0);
    vector<double> sigma_t(number_of_cells * number_of_groups, 1.0);
    vector<double> sigma_s(number_of_cells * number_of_groups * number_of_groups * number_of_scattering_moments, 0.0);
    vector<double> nu_sigma_f(number_of_cells * number_of_groups, 0.0);
    vector<double> chi(number_of_cells * number_of_groups, 0.0);
    vector<string> boundary_conditions(2, "robin"); // dirichlet, marshak, robin
    string problem_type = "adjoint";
    
    for (unsigned i = 0; i < number_of_cells; ++i)
    {
        internal_source[1 + number_of_groups * i] = 1;
        
        for (unsigned gf = 0; gf < number_of_groups; ++gf)
        {
            for (unsigned gt = 0; gt < number_of_groups; ++gt)
            {
                unsigned m = 0;
                unsigned k1 = gf + number_of_groups * (gt + number_of_groups * (i + number_of_cells * m));
                
                m = 1;
                //unsigned k2 = gf + number_of_groups * (gt + number_of_groups * (i + number_of_cells * m));

                if (gf == gt)
                {
                    sigma_s[k1] = 0.4;
                    //sigma_s[k2] = 0.1;
                }
                else if (gf + 1 == gt)
                {
                    sigma_s[k1] = 0.5;
                    //sigma_s[k2] = 0.1;
                }
            }
        }
    }
    
    Transport_Model sp1_transport(number_of_cells,
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
    
    sp1_transport.solve();
    
    sp1_transport.print_scalar_flux();
    
    return 0;
}

int test_spherical()
{
    using namespace std;
    using namespace transport_ns;
    
    unsigned number_of_cells = 3;
    unsigned number_of_groups = 1;
    unsigned number_of_scattering_moments = 1;
    unsigned number_of_ordinates = 8;
    
    unsigned max_iterations = 2;
    double tolerance = 1e-10;
    
    double side_length = 1;//22.017156;
    vector<unsigned> number_of_cells_per_dimension(1, number_of_cells);
    vector<double> internal_source(number_of_cells * number_of_groups, 1.0);
    vector<double> boundary_sources(number_of_groups * number_of_ordinates * 2, 0.0);
    vector<double> sigma_t(number_of_cells * number_of_groups, 0.54628);
    vector<double> sigma_s(number_of_cells * number_of_groups * number_of_groups * number_of_scattering_moments, 0 * 0.464338);
    vector<double> nu_sigma_f(number_of_cells * number_of_groups, 1.70*0.054628);
    vector<double> chi(number_of_cells * number_of_groups, 1.0);
    vector<string> boundary_conditions(2, "reflected");
    boundary_conditions[1] = "reflected";
    string geometry = "spherical";
    
    for (unsigned i = 0; i < number_of_cells; ++i)
    {
        for (unsigned gf = 0; gf < number_of_groups; ++gf)
        {
            for (unsigned gt = 0; gt < number_of_groups; ++gt)
            {
                //unsigned k2 = gf + number_of_groups * (gt + number_of_groups * i);
                
                //sigma_s[k2] = 0.5 / number_of_groups;
            }
        }
    }
    
    Transport_Model sn_transport(number_of_cells,
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
    //sn_transport.solve_eigenvalue();

    sn_transport.print_scalar_flux();
    //sn_transport.print_eigenvalue();
    
    Transport_Model fem_transport(number_of_cells,
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
                                  max_iterations,
                                  tolerance,
                                  geometry);

    fem_transport.solve();
    //fem_transport.solve_eigenvalue();

    fem_transport.print_scalar_flux();
    //fem_transport.print_eigenvalue();
    
    return 0;
}

int test_mc()
{
    using namespace std;
    using namespace mesh_ns;
    using namespace data_ns;
    using namespace monte_carlo_ns;
    
    unsigned number_of_cells = 3;
    unsigned number_of_groups = 1;
    unsigned number_of_scattering_moments = 1;
    unsigned number_of_ordinates = 8;
    
    double side_length = 1;//22.017156;
    vector<double> internal_source(number_of_cells * number_of_groups, 1.0);
    vector<double> boundary_sources(number_of_groups * number_of_ordinates * 2, 0.0);
    vector<double> sigma_t(number_of_cells * number_of_groups, 0.54628);
    vector<double> sigma_s(number_of_cells * number_of_groups * number_of_groups * number_of_scattering_moments, 0 * 0.464338);
    vector<double> nu(number_of_cells * number_of_groups, 1.70);
    vector<double> sigma_f(number_of_cells * number_of_groups, 0.054628);
    vector<double> chi(number_of_cells * number_of_groups, 1.0);
    vector<string> boundary_conditions(2, "reflected");
    boundary_conditions[1] = "reflected";
    
    Mesh mesh(number_of_cells,
              side_length);
    
    Data data(number_of_cells,
              number_of_groups,
              number_of_scattering_moments,
              internal_source,
              boundary_sources,
              sigma_t,
              sigma_s,
              nu,
              sigma_f,
              chi,
              boundary_conditions);

    Monte_Carlo monte_carlo(data,
			    mesh);
    
    return 0;
}

int main(int argc, char *argv[])
{
    //MPI_Init(&argc, &argv);
    
    unsigned checksum = 0;

    // checksum += test_homo();
    // checksum += test_downscatter();
    // checksum += test_adjoint();
    // checksum += test_spherical();
    test_mc();
    
    //MPI_Finalize();
    
    return checksum;
}
