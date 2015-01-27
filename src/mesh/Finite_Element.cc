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

    inline double Finite_Element::
    cell_center_position(unsigned &dimension)
    {
        return cell_center_position_[dimension];
    }
    inline double Finite_Element::
    cell_length(unsigned &dimension)
    {
        return cell_length_[dimension];
    }
}
