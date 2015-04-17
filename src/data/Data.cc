#include "Data.hh"

#include <iostream>
#include <string>
#include <vector>

namespace data_ns
{
    using namespace std;

    Data::
    Data(unsigned &number_of_cells,
         unsigned &number_of_groups,
         unsigned &number_of_scattering_moments,
         vector<double> &internal_source,
         vector<double> &boundary_sources,
         vector<double> &sigma_t,
         vector<double> &sigma_s,
         vector<double> &nu_sigma_f,
         vector<double> &chi,
         vector<string> &boundary_conditions):
        number_of_cells_(number_of_cells),
        number_of_groups_(number_of_groups),
        number_of_scattering_moments_(number_of_scattering_moments),
        internal_source_(internal_source),
        boundary_sources_(boundary_sources),
        sigma_t_(sigma_t),
        sigma_s_(sigma_s),
        sigma_f_(nu_sigma_f),
        chi_(chi),
        boundary_conditions_(boundary_conditions)
    {
        nu_.assign(number_of_groups_ * number_of_cells_, 1);
        
        check();
    }

    Data::
    Data(unsigned &number_of_cells,
         unsigned &number_of_groups,
         unsigned &number_of_scattering_moments,
         vector<double> &internal_source,
         vector<double> &boundary_sources,
         vector<double> &sigma_t,
         vector<double> &sigma_s,
         vector<double> &nu,
         vector<double> &sigma_f,
         vector<double> &chi,
         vector<string> &boundary_conditions):
        number_of_cells_(number_of_cells),
        number_of_groups_(number_of_groups),
        number_of_scattering_moments_(number_of_scattering_moments),
        internal_source_(internal_source),
        boundary_sources_(boundary_sources),
        sigma_t_(sigma_t),
        sigma_s_(sigma_s),
        nu_(nu),
        sigma_f_(sigma_f),
        chi_(chi),
        boundary_conditions_(boundary_conditions)
    {
        check();
    }
    
    int Data::
    check()
    {
        int checksum = 0;
        
        checksum += check_size(internal_source_.size(), number_of_cells_ * number_of_groups_, "internal_source");
        //checksum += check_size(boundary_sources_.size(), number_of_groups_ * 2, "boundary_sources");
        checksum += check_size(sigma_t_.size(), number_of_cells_ * number_of_groups_, "sigma_t");
        checksum += check_size(sigma_s_.size(), number_of_cells_ * number_of_groups_ * number_of_groups_ * number_of_scattering_moments_, "sigma_s");
        checksum += check_size(sigma_f_.size(), number_of_cells_ * number_of_groups_, "nu_sigma_f");
        checksum += check_size(nu_.size(), number_of_cells_ * number_of_groups_, "nu_sigma_f");
        checksum += check_size(chi_.size(), number_of_cells_ * number_of_groups_, "chi");
        checksum += check_size(boundary_conditions_.size(), 2, "boundary_conditions");
        // checksum += check_size(d_.size(), number_of_cells_ * number_of_groups_ * (number_of_scattering_moments_ - 1), "d");

        if (checksum != 0)
        {
            cerr << "data checksum: " << checksum << endl;
        }
        
        return checksum;
    }

    int Data::
    check_size(unsigned vector_size, unsigned expected_size, string vector_name)
    {
        if(vector_size != expected_size)
        {
            cerr << vector_name << " | vector size: " << vector_size << " | expected size: " << expected_size << endl;
            return 1;
        }
        else
        {
            return 0;
        }
    }

} // data_ns
