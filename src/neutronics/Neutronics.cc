#include "Neutronics.hh"

#include <string>
#include <vector>

#include "Data.hh"
#include "Mesh.hh"
#include "Parser.hh"
#include "SPn_Transport.hh"

namespace neutronics_ns
{
    using namespace std;
    using namespace data_ns;
    using namespace mesh_ns;
    using namespace parser_ns;
    using namespace transport_ns;
    
    Neutronics::
    Neutronics(string &input_folder):
        input_folder_(input_folder)
    {
        Parser parser(input_folder_);
        
        parser.parse_data(number_of_dimensions_, "number_of_dimensions");
        parser.parse_data(number_of_cells_, "number_of_cells");
        parser.parse_data(number_of_groups_, "number_of_groups");
        parser.parse_data(number_of_scattering_moments_, "number_of_scattering_moments");
        
        parser.parse_data(side_length_, "side_length");
        parser.parse_data(number_of_cells_per_dimension_, "number_of_cells_per_dimension");
        
        parser.parse_data(internal_source_, "internal_source");
        parser.parse_data(boundary_sources_, "boundary_sources");
        parser.parse_data(sigma_t_, "sigma_t");
        parser.parse_data(sigma_s_, "sigma_s");
        parser.parse_data(nu_sigma_f_, "nu_sigma_f");
        parser.parse_data(chi_, "chi");

        parser.parse_data(boundary_conditions_, "boundary_conditions");

        Neutronics(number_of_dimensions_,
                   number_of_cells_,
                   number_of_groups_,
                   number_of_scattering_moments_,
                   number_of_cells_per_dimension_,
                   side_length_,
                   internal_source_,
                   boundary_sources_,
                   sigma_t_,
                   sigma_s_,
                   nu_sigma_f_,
                   chi_,
                   boundary_conditions_);
    }

    Neutronics::
    Neutronics(unsigned &number_of_dimensions,
               unsigned &number_of_cells,
               unsigned &number_of_groups,
               unsigned &number_of_scattering_moments,
               vector<unsigned> &number_of_cells_per_dimension,
               vector<double> &side_length,
               vector<double> &internal_source,
               vector<double> &boundary_sources,
               vector<double> &sigma_t,
               vector<double> &sigma_s,
               vector<double> &nu_sigma_f,
               vector<double> &chi,
               vector<string> &boundary_conditions)
    {
        mesh_ = new Mesh(number_of_dimensions,
                         number_of_cells_per_dimension,
                         side_length);
        
        data_ = new Data(number_of_cells,
                         number_of_groups,
                         number_of_scattering_moments,
                         internal_source,
                         boundary_sources,
                         sigma_t,
                         sigma_s,
                         nu_sigma_f,
                         chi,
                         boundary_conditions);

        unsigned number_of_moments = (number_of_scattering_moments - 1) / 2;
        
        transport_ = new SPn_Transport(number_of_moments,
                                       *data_,
                                       *mesh_);
    }
    
    Neutronics::
    ~Neutronics()
    {
        delete data_;
        delete mesh_;
        delete transport_;
    }

    void Neutronics::
    solve()
    {
        transport_->solve();
    }
}
