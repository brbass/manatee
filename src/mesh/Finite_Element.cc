#include "Finite_Element.hh"

#include <vector>

#include <iostream>
using namespace std;

namespace mesh_ns
{
    Finite_Element::
    Finite_Element(vector<double> cell_center_position,
                   vector<double> cell_length):
        cell_center_position_(cell_center_position),
        cell_length_(cell_length)
    {
    }
}
