#ifndef Data_hh
#define Data_hh

#include <string>
#include <vector>

namespace data_ns
{
    using std::string;
    using std::vector;
    
    class Data
    {
    private:
        
        unsigned number_of_cells_;
        unsigned number_of_groups_;
        unsigned number_of_scattering_moments_;
        
        vector<double> &internal_source_;
        vector<double> &boundary_sources_;
        vector<double> &sigma_t_;
        vector<double> &sigma_s_;
        vector<double> &nu_sigma_f_;
        vector<double> &chi_;
        vector<double> d_;
        
        vector<string> &boundary_conditions_;

        int check_size(unsigned vector_size, unsigned expected_size, string vector_name);

    public:
        
        Data(unsigned &number_of_cells,
             unsigned &number_of_groups,
             unsigned &number_of_scattering_moments,
             vector<double> &internal_source,
             vector<double> &boundary_sources,
             vector<double> &sigma_t,
             vector<double> &sigma_s,
             vector<double> &nu_sigma_f,
             vector<double> &chi,
             vector<string> &boundary_conditions);

        //void compute_d();
        int check();
        
        inline unsigned number_of_groups()
        {
            return number_of_groups_;
        }
    
        inline unsigned number_of_scattering_moments()
        {
            return number_of_scattering_moments_;
        }
    
        inline double internal_source(unsigned cell, unsigned group)
        {
            return internal_source_[group + number_of_groups_ * cell];
        }
    
        inline double boundary_source(unsigned boundary, unsigned group, unsigned ordinate = 0)
        {
            return boundary_sources_[group + number_of_groups_ * (boundary + 2 * ordinate)];
        }
        
        inline double sigma_t(unsigned cell, unsigned group)
        {
            return sigma_t_[group + number_of_groups_ * cell];
        }
    
        inline double sigma_s(unsigned cell, unsigned from_group, unsigned to_group, unsigned moment = 0)
        {
            return sigma_s_[from_group + number_of_groups_ * (to_group + number_of_groups_ * (cell + number_of_cells_ * moment))];
        }
    
        inline double nu_sigma_f(unsigned cell, unsigned from_group)
        {
            return nu_sigma_f_[from_group + number_of_groups_ * cell];
        }
            
        inline double chi(unsigned cell, unsigned to_group)
        {
            return chi_[to_group * number_of_groups_ * cell];
        }
    
        inline double d(unsigned cell, unsigned from_group, unsigned moment)
        {
            return d_[from_group + number_of_groups_ * (cell + number_of_cells_ * moment)];
        }

        inline string boundary_condition(unsigned boundary)
        {
            return boundary_conditions_[boundary];
        }

    }; // Data
} // data_ns

#endif // Data_hh
