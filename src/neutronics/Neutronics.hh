#ifndef Neutronics_hh
#define Neutronics_hh

#include <memory>
#include <string>
#include <vector>

#include "Data.hh"
#include "Mesh.hh"
#include "SP1_Transport.hh"
#include "SPn_Transport.hh"

namespace neutronics_ns
{
    using std::string;
    using std::unique_ptr;
    using std::vector;
    using namespace data_ns;
    using namespace mesh_ns;
    using namespace transport_ns;
    
    class Neutronics
    {
        
    private:

        unique_ptr<Data> data_;
        unique_ptr<Mesh> mesh_;
        unique_ptr<SP1_Transport> sp1_transport_;
        unique_ptr<SPn_Transport> spn_transport_;
        
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
        
        string input_folder_;
        
    public:

        Neutronics(string &input_folder);

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
                   vector<string> &boundary_conditions);

        void solve();

        void print_scalar_flux()
        {
            //spn_transport_->print_scalar_flux();
            sp1_transport_->print_scalar_flux();
        }
    };
}

#endif
