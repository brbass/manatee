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

        void compute_d();
        
        inline unsigned number_of_groups();
        inline unsigned number_of_scattering_moments();
        
        inline double internal_source(unsigned &cell, unsigned &group);
        inline double boundary_source(unsigned &ordinate, unsigned &group);
        inline double sigma_t(unsigned &cell, unsigned &group);
        inline double sigma_s(unsigned &cell, unsigned &from_group, unsigned &to_group, unsigned &moment);
        inline double nu_sigma_f(unsigned &cell, unsigned &from_group);
        inline double chi(unsigned &cell, unsigned &to_group);
        inline double d(unsigned &cell, unsigned &from_group);
        
        inline string boundary_condition(unsigned &boundary);
        
    }; // Data
} // data_ns

#endif // Data_hh
