#ifndef Mesh_hh
#define Mesh_hh

#include <vector>

#include "Finite_Element.hh"

namespace mesh_ns
{
    using std::vector;

    class Mesh
    {

    private:
    
        unsigned &number_of_dimensions_;
        unsigned number_of_cells_;
        unsigned number_of_nodes_ = 2;
        
        vector<Finite_Element> elements_;
        
    public:
    
        Mesh(unsigned &number_of_dimensions,
             vector<unsigned> &number_of_cells_per_dimension,
             vector<double> &side_length);

        int check();
        
        // Mesh(unsigned &number_of_dimensions,
        //      unsigned number_of_cells,
        //      vector<double> &cell_boundaries);

        inline unsigned number_of_dimensions()
        {
            return number_of_dimensions_;
        }
        
        inline unsigned number_of_cells()
        {
            return number_of_cells_;
        }
        
        inline unsigned number_of_nodes()
        {
            return number_of_nodes_;
        }

        inline double cell_length(unsigned cell)
        {
            return elements_[cell].cell_length(0);
        }

        inline double stiffness(unsigned cell, unsigned a, unsigned b)
        {
            return elements_[cell].stiffness(a, b);
        }

        inline double stiffness_moment(unsigned cell, unsigned a, unsigned b)
        {
            return elements_[cell].stiffness_moment(a, b);
        }
    };
}

#endif
