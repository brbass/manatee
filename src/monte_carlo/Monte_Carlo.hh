#ifndef Monte_Carlo_hh
#define Monte_Carlo_hh

#include <random>

#include "Data.hh"
#include "Mesh.hh"

namespace monte_carlo_ns
{
    using namespace data_ns;
    using namespace mesh_ns;
 
    class Monte_Carlo
    {   
    public:
    
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
                vector<string> boundary_conditions);
        
    private:

        enum reaction
        {
            ABSORPTION,
            SCATTERING,
            FISSION,
            LEAKAGE
        };
        
        Data data_;
        Mesh mesh_;
        
        std::mt19937 generator_;
        std::uniform_real_distribution<double> distribution_;
        
        inline double get_rand()
        {
            return distribution_(generator_);
        }

        double sample_source();
        
        double sample_reaction();
        
        
    };
}

#endif
