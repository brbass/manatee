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
        
        inline double cell_center_position(unsigned dimension)
        {
            return cell_center_position_[dimension];
        }
        
        inline double cell_length()
        {
            return cell_length_[0];
        }

        inline double stiffness(unsigned a, unsigned b)
        {
            if (a == b)
            {
                return 2.0 / 3.0;
            }
            else
            {
                return 1.0 / 3.0;
            }
        }
    
        inline double stiffness_moment(unsigned a, unsigned b)
        {
            if (a == b)
            {
                return 1.0 / 2.0;
            }
            else
            {
                return -1.0 / 2.0;
            }
        }   
    };
}

#endif
