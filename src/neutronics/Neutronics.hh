#ifndef Neutronics_hh
#define Neutronics_hh

#include <string>
#include <vector>

#include "../data/Data.hh"
#include "../mesh/Mesh.hh"
#include "../parser/Parser.hh"

namespace neutronics_ns
{
    using std::string;
    using std::vector;
    using namespace data_ns;
    using namespace mesh_ns;
    using namespace parser_ns;

    
    class Neutronics
    {
        
    private:

        Data *data;
        Mesh *mesh;
        Parser *parser;

        unsigned number_of_dimensions_;
        unsigned number_of_cells_;
        unsigned number_of_groups_;
        unsigned number_of_scattering_moments_;

        vector<unsigned> number_of_cells_per_dimension_;

        vector<double> side_length_;
        vector<double> internal_source_;
        vector<double> boundary_sources_;
        vector<double> sigma_t_;
        vector<double> sigma_s_;
        vector<double> nu_sigma_f_;
        vector<double> chi_;

        vector<string> boundary_conditions_;
        
        string &input_folder_;
        
    public:

        Neutronics(string &input_folder);

        ~Neutronics();
        
    };
}

#endif
