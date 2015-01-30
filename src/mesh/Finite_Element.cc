#include "Finite_Element.hh"

#include <vector>

namespace mesh_ns
{
    Finite_Element::
    Finite_Element(vector<double> cell_center_position,
                   vector<double> cell_length):
        cell_center_position_(cell_center_position),
        cell_length_(cell_length)
    {
    }

    double Finite_Element::
    stiffness(unsigned a, unsigned b)
    {
        if (a == b)
        {
            return 2 / 3;
        }
        else
        {
            return 1 / 3;
        }
    }
    
    double Finite_Element::
    stiffness_moment(unsigned a, unsigned b)
    {
        if (a == b)
        {
            return 1 / 2;
        }
        else
        {
            return -1 / 2;
        }
    }   
}
