#include "Finite_Element.hh"

namespace mesh_ns
{
    Finite_Element::
    Finite_Element(double cell_center_position,
                   double cell_length):
        cell_center_position_(cell_center_position),
        cell_length_(cell_length)
    {
    }
}
