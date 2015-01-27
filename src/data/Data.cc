#include "Data.hh"

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
    }

    void Data::
    compute_d()
    {
    }
    
    inline unsigned Data::
    number_of_groups()
    {
        return number_of_groups_;
    }
    
    inline unsigned Data::
    number_of_scattering_moments()
    {
        return number_of_scattering_moments_;
    }
    
    inline double Data::
    internal_source(unsigned &cell, unsigned &group)
    {
        return internal_source_[group + number_of_groups_ * cell];
    }
    
    inline double Data::
    boundary_source(unsigned &ordinate, unsigned &group)
    {
        return boundary_sources_[group + number_of_groups_ * ordinate];
    }
    
    inline double Data::
    sigma_t(unsigned &cell, unsigned &group)
    {
        return sigma_t_[group + number_of_groups_ * cell];
    }
    
    inline double Data::
    sigma_s(unsigned &cell, unsigned &from_group, unsigned &to_group, unsigned &moment)
    {
        return sigma_s_[from_group + number_of_groups_ * (to_group + number_of_groups_ * (cell + number_of_cells_ * moment))];
    }
    
    inline double Data::
    nu_sigma_f(unsigned &cell, unsigned &from_group)
    {
        return nu_sigma_f_[from_group + number_of_groups_ * cell];
    }
            
    inline double Data::
    chi(unsigned &cell, unsigned &to_group)
    {
        return chi_[to_group * number_of_groups_ * cell];
    }
    
    inline double Data::
    d(unsigned &cell, unsigned &from_group)
    {
        return d_[from_group + number_of_groups_ * cell];
    }

    inline string Data::
    boundary_condition(unsigned &boundary)
    {
        return boundary_conditions_[boundary];
    }
    
} // data_ns
