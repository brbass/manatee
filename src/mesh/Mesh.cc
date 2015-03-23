#include "Mesh.hh"

#include <vector>
#include <iostream>

#include "Finite_Element.hh"

namespace mesh_ns
{
    using namespace std;
    
    Mesh::
    Mesh(unsigned &number_of_cells,
         double &side_length):
        number_of_cells_(number_of_cells)
    {
        elements_.resize(0);
        
        double average_cell_length = side_length / number_of_cells_;
        double current_position = average_cell_length / 2;
        
        for (unsigned i = 0; i < number_of_cells_; ++i)
        {
            elements_.emplace_back(current_position,
                                   average_cell_length;
                
            current_position += average_cell_length;
        }
            
        check();
    }

    int Mesh::
    check()
    {
        Int checksum = 0;
        
        checksum += check_size(elements_.size(), number_of_cells_, "number_of_cells");
        
        if (checksum != 0)
        {
            cout << "mesh checksum: " << checksum << endl;
        }

        return checksum;
    }

    int Mesh::
    check_size(unsigned vector_size, unsigned expected_size, string vector_name)
    {
        if(vector_size != expected_size)
        {
            cout << vector_name << " | vector size: " << vector_size << " | expected size: " << expected_size << endl;
            return 1;
        }
        else
        {
            return 0;
        }
    }

}
