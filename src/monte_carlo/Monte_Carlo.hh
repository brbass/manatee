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

        Data &data_;
        Mesh &mesh_;

        std::mt19937 generator_;
        std::uniform_real_distribution<double> distribution_;

        inline double get_rand()
        {
	     return distribution_(generator_);
        }
    };
}
