#ifndef Mesh_hh
#define Mesh_hh

#include <iostream>
using namespace std;
#include <string>
#include <vector>

#include "Finite_Element.hh"

namespace mesh_ns
{
    using std::string;
    using std::vector;

    class Mesh
    {

    private:
    
        unsigned number_of_cells_;
        unsigned number_of_nodes_ = 2;
        unsigned side_length_;
        string geometry_;
        
        vector<Finite_Element> elements_;
        
        int check_size(unsigned vector_size, unsigned expected_size, string vector_name);
        
    public:
    
        Mesh(unsigned &number_of_cells,
             double &side_length,
             string geometry = "slab");

        int check();
        
        // Mesh(unsigned &number_of_dimensions,
        //      unsigned number_of_cells,
        //      vector<double> &cell_boundaries);

        inline unsigned number_of_cells()
        {
            return number_of_cells_;
        }
        
        inline unsigned number_of_nodes()
        {
            return number_of_nodes_;
        }
        
        inline double cell_center_position(unsigned cell)
        {
            return elements_[cell].cell_center_position();
        }

        inline double cell_edge_position(unsigned cell, unsigned edge)
        {
            return elements_[cell].cell_edge_position(edge);
        }
        
        inline double cell_length(unsigned cell)
        {
            return elements_[cell].cell_length();
        }
        
        inline double stiffness(unsigned cell, unsigned kb, unsigned kw, unsigned kdb = 0, unsigned kdw = 0)
        {
            return elements_[cell].stiffness(kb, kw, kdb, kdw);
        }
        
        inline bool geometry(string test)
        {
            return geometry_.compare(test) == 0;
        }

        inline double side_length()
        {
            return side_length_;
        }
    };
}

#endif
