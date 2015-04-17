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
    
        Monte_Carlo(Data &data,
                    Mesh &mesh);
        
    private:

        enum reaction
        {
            ABSORPTION,
            SCATTERING,
            FISSION,
            LEAKAGE
        };
        
        Data &data_;
        Mesh &mesh_;
        
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
