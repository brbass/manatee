#ifndef Transport_Model_hh
#define Transport_Model_hh

#include <memory>
#include <string>
#include <vector>

#include "Data.hh"
#include "FEM_Transport.hh"
#include "Mesh.hh"
#include "Ordinates.hh"
#include "Sn_Transport.hh"
#include "SP1_Transport.hh"
#include "SPn_Transport.hh"

namespace transport_ns
{
    using std::string;
    using std::unique_ptr;
    using std::vector;
    using namespace data_ns;
    using namespace mesh_ns;
    
    class Transport_Model
    {
        
    private:

        unique_ptr<Data> data_;
        unique_ptr<FEM_Transport> fem_transport_;
        unique_ptr<Mesh> mesh_;
        unique_ptr<Ordinates> ordinates_;
        unique_ptr<SP1_Transport> sp1_transport_;
        unique_ptr<SPn_Transport> spn_transport_;
        unique_ptr<Sn_Transport> sn_transport_;
        
        unsigned number_of_dimensions_;
        unsigned number_of_cells_;
        unsigned number_of_groups_;
        unsigned number_of_scattering_moments_;

        double side_length_;

        string geometry_;
        string transport_type_;
        string problem_type_;
        
        vector<double> internal_source_;
        vector<double> boundary_sources_;
        vector<double> sigma_t_;
        vector<double> sigma_s_;
        vector<double> nu_sigma_f_;
        vector<double> chi_;

        vector<string> boundary_conditions_;
        
        string input_folder_;
        
    public:

        Transport_Model(string &input_folder);

        Transport_Model(unsigned &number_of_cells,
                        unsigned &number_of_groups,
                        unsigned &number_of_scattering_moments,
                        double &side_length,
                        vector<double> &internal_source,
                        vector<double> &boundary_sources,
                        vector<double> &sigma_t,
                        vector<double> &sigma_s,
                        vector<double> &nu_sigma_f,
                        vector<double> &chi,
                        vector<string> &boundary_conditions,
                        string problem_type = "forward");

        Transport_Model(unsigned &number_of_cells,
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
                        vector<string> &boundary_conditions,
                        string geometry = "slab");

        Transport_Model(unsigned &number_of_cells,
                        unsigned &number_of_groups,
                        unsigned &number_of_scattering_moments,
                        double &side_length,
                        vector<double> &internal_source,
                        vector<double> &boundary_sources,
                        vector<double> &sigma_t,
                        vector<double> &sigma_s,
                        vector<double> &nu_sigma_f,
                        vector<double> &chi,
                        vector<string> &boundary_conditions,
                        unsigned &max_iterations,
                        double &tolerance,
                        string geometry = "spherical");
        
        void solve();
        void solve_eigenvalue();
        
        void print_scalar_flux()
        {
            if (transport_type_.compare("sn_transport") == 0)
            {
                sn_transport_->print_scalar_flux();
            }
            else if (transport_type_.compare("sp1_transport") == 0)
            {
                sp1_transport_->print_scalar_flux();
            }
            else if (transport_type_.compare("fem_transport") == 0)
            {
                fem_transport_->print_scalar_flux();
            }
            else
            {
                cerr << "transport type does not exist" << endl;
            }
        }

        void print_eigenvalue()
        {
            if (transport_type_.compare("sn_transport") == 0)
            {
                sn_transport_->print_eigenvalue();
            }
            // else if (transport_type_.compare("sp1_transport") == 0)
            // {
            //     sp1_transport_->print_scalar_flux();
            // }
            else if (transport_type_.compare("fem_transport") == 0)
            {
                fem_transport_->print_eigenvalue();
            }
            else
            {
                cerr << "transport type does not exist" << endl;
            }
        }

    };
}

#endif
