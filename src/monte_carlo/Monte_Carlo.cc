#include "Monte_Carlo.hh"

#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "Data.hh"
#include "Mesh.hh"

namespace monte_carlo_ns
{
    using namespace std;
    using namespace data_ns;
    using namespace mesh_ns;
    
    Monte_Carlo::
    Monte_Carlo(unsigned &number_of_cells,
                unsigned &number_of_groups,
                unsigned &number_of_scattering_moments,
                double &side_length,
                vector<double> &internal_source,
                vector<double> &boundary_sources,
                vector<double> &sigma_t,
                vector<double> &sigma_s,
                vector<double> &nu,
                vector<double> &sigma_f,
                vector<double> &chi,
                vector<string> boundary_conditions):
        data_(number_of_cells,
              number_of_groups,
              number_of_scattering_moments,
              internal_source,
              boundary_sources,
              sigma_t,
              sigma_s,
              nu,
              sigma_f,
              chi,
              boundary_conditions),
        mesh_(number_of_cells,
              side_length),
	distribution_(0.0, 1.0)
    {
        data_.check();
        mesh_.check();
    }
}
