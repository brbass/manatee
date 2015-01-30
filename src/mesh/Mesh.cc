#include "Mesh.hh"

#include <vector>
#include <iostream>

#include "Finite_Element.hh"

namespace mesh_ns
{
    using namespace std;
    
    Mesh::
    Mesh(unsigned &number_of_dimensions,
         vector<unsigned> &number_of_cells_per_dimension,
         vector<double> &side_length):
        number_of_dimensions_(number_of_dimensions)
    {
        elements_.resize(0);
        
        if (number_of_dimensions_ == 1)
        {
            number_of_cells_ = number_of_cells_per_dimension[0];

            double average_cell_length = side_length[0] / number_of_cells_;
            double current_position = average_cell_length / 2;
            
            for (unsigned i = 0; i < number_of_cells_; ++i)
            {
                vector<double> cell_center_position(current_position, 1);
                vector<double> cell_length(average_cell_length, 1);
                
                Finite_Element element(cell_center_position,
                                       cell_length);

                elements_.push_back(element);
                
                current_position += average_cell_length;
            }
        }
        else
        {
            cout << "Only one dimension is supported" << endl;
        }
    }

    int Mesh::
    check()
    {
    }
    
    inline unsigned Mesh::
    number_of_dimensions()
    {
        return number_of_dimensions_;
    }
    inline unsigned Mesh::
    number_of_cells()
    {
        return number_of_cells_;
    }
    inline unsigned Mesh::
    number_of_nodes()
    {
        return number_of_nodes_;
    }

    inline double Mesh::
    cell_length(unsigned &cell)
    {
        return elements_[cell].cell_length();
    }

    inline double Mesh::
    stiffness(unsigned &cell, unsigned a, unsigned b)
    {
        return elements_[cell].stiffness(a, b);
    }

    inline double Mesh::
    stiffness_moment(unsigned &cell, unsigned a, unsigned b)
    {
        return elements_[cell].stiffness_moment(a, b);
    }
}
