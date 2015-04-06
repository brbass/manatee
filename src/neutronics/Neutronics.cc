#include "Neutronics.hh"

#include <iostring>
#include <string>
#include <vector>

#include "Data.hh"
#include "Mesh.hh"
#include "Ordinates.hh"
#include "Parser.hh"
#include "Sn_Transport.hh"
#include "SP1_Transport.hh"
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
        
        parser.parse_data(internal_source_, "internal_source");
        parser.parse_data(boundary_sources_, "boundary_sources");
        parser.parse_data(sigma_t_, "sigma_t");
        parser.parse_data(sigma_s_, "sigma_s");
        parser.parse_data(nu_sigma_f_, "nu_sigma_f");
        parser.parse_data(chi_, "chi");
        
        parser.parse_data(boundary_conditions_, "boundary_conditions");
        
        Neutronics(number_of_cells_,
                   number_of_groups_,
                   number_of_scattering_moments_,
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
    Neutronics(unsigned &number_of_cells,
               unsigned &number_of_groups,
               unsigned &number_of_scattering_moments,
               double &side_length,
               vector<double> &internal_source,
               vector<double> &boundary_sources,
               vector<double> &sigma_t,
               vector<double> &sigma_s,
               vector<double> &nu_sigma_f,
               vector<double> &chi,
               vector<string> &boundary_conditions)
    {
        mesh_ = unique_ptr<Mesh> (new Mesh(number_of_cells,
                                           side_length));
        
        data_ = unique_ptr<Data> (new Data(number_of_cells,
                                           number_of_groups,
                                           number_of_scattering_moments,
                                           internal_source,
                                           boundary_sources,
                                           sigma_t,
                                           sigma_s,
                                           nu_sigma_f,
                                           chi,
                                           boundary_conditions));
        
        // unsigned number_of_moments = number_of_scattering_moments / 2;
        
        // spn_transport_ = unique_ptr<SPn_Transport> (new SPn_Transport(number_of_moments,
        //                                                              *data_,
        //                                                              *mesh_));

        sp1_transport_ = unique_ptr<SP1_Transport> (new SP1_Transport(*data_,
                                                                      *mesh_));

        transport_type_ = "sn_transport";
    }

    Neutronics::
    Neutronics(unsigned &number_of_cells,
               unsigned &number_of_groups,
               unsigned &number_of_scattering_moments,
               unsigned &number_of_ordinates,
               double &side_length,
               vector<double> &internal_source,
               vector<double> &boundary_sources,
               vector<double> &sigma_t,
               vector<double> &sigma_s,
               vector<double> &nu_sigma_f,
               vector<double> &chi,
               vector<string> &boundary_conditions)
    {
        mesh_ = unique_ptr<Mesh> (new Mesh(number_of_cells,
                                           side_length));
        
        data_ = unique_ptr<Data> (new Data(number_of_cells,
                                           number_of_groups,
                                           number_of_scattering_moments,
                                           internal_source,
                                           boundary_sources,
                                           sigma_t,
                                           sigma_s,
                                           nu_sigma_f,
                                           chi,
                                           boundary_conditions));
        
        ordinates_ = unique_ptr<Ordinates> (new Ordinates(number_of_ordinates));

        sn_transport_ = unique_ptr<Sn_Transport> (new Sn_Transport(*data_,
                                                                   *mesh_,
                                                                   *ordinates_));

        transport_type_ = "sn_transport";
    }
    
    void Neutronics::
    solve()
    {
        // spn_transport_->solve();
        if (transport_type_.compare("sn_transport") == 0)
        {
            sn_transport_->solve();
        }
        else if (transport_type_.compare("sp1_transport") == 0)
        {
            sp1_transport_->solve();
        }
        else
        {
            cerr << "transport type does not exist" << endl;
        }
    }
}
