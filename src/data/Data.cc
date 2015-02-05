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
        nu_sigma_f_(nu_sigma_f),
        chi_(chi),
        boundary_conditions_(boundary_conditions)
    {
        compute_d();
        int checksum = check();
    }

    void Data::
    compute_d()
    {
        d_.resize(number_of_cells_ * number_of_groups_ * (number_of_scattering_moments_ - 1));

        for (unsigned i = 0; i < number_of_cells_; ++i)
        {
            for (unsigned g = 0; g < number_of_groups_; ++g)
            {
                for (unsigned m = 0; m < number_of_scattering_moments_ - 1; ++m)
                {
                    double sum = 0;

                    for (unsigned g2 = 0; g2 < number_of_groups_; ++g2)
                    {
                        unsigned k = g + number_of_groups_ * (g2 + number_of_groups_ * (i + number_of_cells_ * (m + 1)));
                        
                        sum += sigma_s_[k];
                    }
                    
                    unsigned k1 = g + number_of_groups_ * (i + number_of_cells_ * m);
                    unsigned k2 = g + number_of_groups_ * i;

                    d_[k1] = 1.0 * ((m+1)*(m+1))/((2*m+1)*(2*m+3)) * (sigma_t_[k2] - sum);
                }
            }
        }
    }

    int Data::
    check()
    {
        int checksum = 0;
        
        checksum += check_size(internal_source_.size(), number_of_cells_ * number_of_groups_, "internal_source");
        checksum += check_size(boundary_sources_.size(), number_of_scattering_moments_ * 2, "boundary_sources");
        checksum += check_size(sigma_t_.size(), number_of_cells_ * number_of_groups_, "sigma_t");
        checksum += check_size(sigma_s_.size(), number_of_cells_ * number_of_groups_ * number_of_groups_ * number_of_scattering_moments_, "sigma_s");
        checksum += check_size(nu_sigma_f_.size(), number_of_cells_ * number_of_groups_, "nu_sigma_f");
        checksum += check_size(chi_.size(), number_of_cells_ * number_of_groups_, "chi");
        checksum += check_size(boundary_conditions_.size(), 2, "boundary_conditions");
        checksum += check_size(d_.size(), number_of_cells_ * number_of_groups_ * (number_of_scattering_moments_ - 1), "d");

        if (checksum != 0)
        {
            cout << "data checksum: " << checksum << endl;
        }
        
        return checksum;
    }

    int Data::
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

} // data_ns
