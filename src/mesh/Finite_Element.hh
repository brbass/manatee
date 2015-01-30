#ifndef Finite_Element_hh
#define Finite_Element_hh

#include <vector>

namespace mesh_ns
{
    using std::vector;

    class Finite_Element
    {
        
    private:
        
        vector<double> cell_center_position_;
        vector<double> cell_length_;
        
    public:

        Finite_Element(){}
        
        Finite_Element(vector<double> cell_center_position,
                       vector<double> cell_length_);
        
        inline double cell_center_position(unsigned &dimension);
        inline double cell_length();//unsigned &dimension);
        double stiffness(unsigned a, unsigned b);
        double stiffness_moment(unsigned a, unsigned b);
    };
}

#endif
