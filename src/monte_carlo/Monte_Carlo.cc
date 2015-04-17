#include "Monte_Carlo.hh"

#include <iostream>
#include <random>
#include <string>

#include "Data.hh"
#include "Mesh.hh"

namespace monte_carlo_ns
{
    using namespace std;
    using namespace data_ns;
    using namespace mesh_ns;
    
    Monte_Carlo::
    Monte_Carlo(Data &data,
                Mesh &mesh)
        :
        data_(data),
        mesh_(mesh),
	distribution_(0.0, 1.0)
    {
        data_.check();
        mesh_.check();
    }
}
