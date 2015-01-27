#include "Neutronics.hh"

#include <string>
#include <vector>

#include "../data/Data.hh"
#include "../mesh/Mesh.hh"
#include "../parser/Parser.hh"

namespace neutronics_ns
{
    using namespace std;
    using namespace data_ns;
    using namespace mesh_ns;
    using namespace parser_ns;
    
    Neutronics::
    Neutronics(string &input_folder):
        input_folder_(input_folder)
    {
        parser = new Parser(input_folder_);

        parser->parse_data(number_of_dimensions_, "number_of_dimensions");
        parser->parse_data(number_of_cells_, "number_of_cells");
        parser->parse_data(number_of_groups_, "number_of_groups");
        parser->parse_data(number_of_scattering_moments_, "number_of_scattering_moments");
        
        parser->parse_data(side_length_, "side_length");
        parser->parse_data(number_of_cells_per_dimension_, "number_of_cells_per_dimension");
        
        parser->parse_data(internal_source_, "internal_source");
        parser->parse_data(boundary_sources_, "boundary_sources");
        parser->parse_data(sigma_t_, "sigma_t");
        parser->parse_data(sigma_s_, "sigma_s");
        parser->parse_data(nu_sigma_f_, "nu_sigma_f");
        parser->parse_data(chi_, "chi");

        parser->parse_data(boundary_conditions_, "boundary_conditions");
        
        mesh = new Mesh(number_of_dimensions_,
                        number_of_cells_per_dimension_,
                        side_length_);
        
        data = new Data(number_of_cells_,
                        number_of_groups_,
                        number_of_scattering_moments_,
                        internal_source_,
                        boundary_sources_,
                        sigma_t_,
                        sigma_s_,
                        nu_sigma_f_,
                        chi_,
                        boundary_conditions_);
    }

    Neutronics::
    ~Neutronics()
    {
        delete data;
        delete mesh;
        delete parser;
    }
}
