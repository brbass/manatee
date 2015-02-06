#include <iostream>
#include <string>

#include "mpi.h"

#include "Neutronics.hh"

int test_neutronics()
{
    using namespace std;
    using namespace neutronics_ns;
    
    unsigned checksum = 0;
    
    unsigned number_of_dimensions = 1;
    unsigned number_of_cells = 10;
    unsigned number_of_groups = 3;
    unsigned number_of_scattering_moments = 3;
    
    vector<double> side_length(1, 2.0);
    vector<unsigned> number_of_cells_per_dimension(1, number_of_cells);
    vector<double> internal_source(number_of_cells * number_of_groups, 1.0);
    vector<double> boundary_sources(number_of_scattering_moments * 2, 0.0);
    vector<double> sigma_t(number_of_cells * number_of_groups, 1.0);
    vector<double> sigma_s(number_of_cells * number_of_groups * number_of_groups * number_of_scattering_moments, 0.0);
    vector<double> nu_sigma_f(number_of_cells * number_of_groups, 0.0);
    vector<double> chi(number_of_cells * number_of_groups, 0.0);
    vector<string> boundary_conditions(2, "vacuum");

    for (unsigned i = 0; i < number_of_cells; ++i)
    {
        for (unsigned g1 = 0; g1 < number_of_groups; ++g1)
        {
            for (unsigned g2 = 0; g2 < number_of_groups; ++g2)
            {
                unsigned m = 0;
                unsigned k = g1 + number_of_groups * (g2 + number_of_groups * (i + number_of_cells * m));
                
                sigma_s[k] = 0.2 / number_of_groups;

                m = 1;
                k = g1 + number_of_groups * (g2 + number_of_groups * (i + number_of_cells * m));

                sigma_s[k] = 0.1 / number_of_groups;

                m = 2;
                k = g1 + number_of_groups * (g2 + number_of_groups * (i + number_of_cells * m));

                sigma_s[k] = 0.05 / number_of_groups;
            }
        }
    }
    
    Neutronics neutronics(number_of_dimensions,
                          number_of_cells,
                          number_of_groups,
                          number_of_scattering_moments,
                          number_of_cells_per_dimension,
                          side_length,
                          internal_source,
                          boundary_sources,
                          sigma_t,
                          sigma_s,
                          nu_sigma_f,
                          chi,
                          boundary_conditions);
    
    neutronics.solve();

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