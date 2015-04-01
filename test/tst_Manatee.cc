#include <iostream>
#include <string>

#include "mpi.h"

#include "Neutronics.hh"

int test_neutronics()
{
    using namespace std;
    using namespace neutronics_ns;
    
    unsigned checksum = 0;
    
    unsigned number_of_cells = 10;
    unsigned number_of_groups = 3;
    unsigned number_of_scattering_moments = 2;
    
    double side_length = 1;
    vector<unsigned> number_of_cells_per_dimension(1, number_of_cells);
    vector<double> internal_source(number_of_cells * number_of_groups, 1.0);
    vector<double> boundary_sources(number_of_groups * 2, 1.0);
    vector<double> sigma_t(number_of_cells * number_of_groups, 0.5);
    vector<double> sigma_s(number_of_cells * number_of_groups * number_of_groups * number_of_scattering_moments, 0);
    vector<double> nu_sigma_f(number_of_cells * number_of_groups, 0.0);
    vector<double> chi(number_of_cells * number_of_groups, 0.0);
    vector<string> boundary_conditions(2, "robin"); // dirichlet, marshak, robin
    
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
    
    // for (unsigned i = 0; i < number_of_cells; ++i)
    // {
    //     for (unsigned gf = 0; gf < number_of_groups; ++gf)
    //     {
    //         if (gf == 0)
    //         {
    //             internal_source[gf + number_of_groups * i] = 1.0 * (i + number_of_cells) / number_of_cells;
    //         }

    //         unsigned k1 = gf + number_of_groups * i;
            
    //         sigma_t[k1] *= 0.2 * (i + number_of_cells) / number_of_cells;
            
    //         for (unsigned gt = 0; gt < number_of_groups; ++gt)
    //         {
    //             for (unsigned m = 0; m < number_of_scattering_moments; ++m)
    //             {
    //                 unsigned k2 = gf + number_of_groups * (gt + number_of_groups * (i + number_of_cells * m));
    //                 sigma_s[k2] = 0.2 * sigma_t[k1] / number_of_groups / (m + 1);
    //                 sigma_s[k2] = 0.2 * (gf + 1) * (2 * gf + 1) / number_of_groups / number_of_groups / (m + 1) / number_of_groups;
    //                 sigma_s[k2] *= 1.0 * i / number_of_cells;
                    
    //                 sigma_s[k2] = 2.0 / number_of_groups;
                    
    //                 if (gt == gf + 1)
    //                 {
    //                     sigma_s[k2] = 0.5;
    //                 }
    //             }
    //         }
    //     }
    // }
    
    Neutronics neutronics(number_of_cells,
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
    
    neutronics.solve();
    
    neutronics.print_scalar_flux();
    
    return checksum;
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    
    unsigned checksum = 0;
    
    checksum += test_neutronics();
    
    MPI_Finalize();
    
    return checksum;
}
