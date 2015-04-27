#ifndef Monte_Carlo_hh
#define Monte_Carlo_hh

#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#include "Bank.hh"
#include "Data.hh"
#include "Estimators.hh"
#include "Mesh.hh"

namespace monte_carlo_ns
{
    using std::vector;
    using namespace data_ns;
    using namespace mesh_ns;
 
    class Monte_Carlo
    {   
    private:
        
        struct Weight_Window
        {
            double upper = 2;
            double ideal = 1;
            double lower = 0.5;
        };
        
        bool implicit_capture_;
        bool init_with_adjoint_ = false;
        unsigned number_of_histories_;

        double max_splitting_ = 1.0;
        double runtime_ = 0.0;
        
        double total_source_;
        // double total_birth_weight_ = 0.0;
        double birth_weight_ = 1.0;
        double upper_over_lower_ = 4.0;
        double ideal_over_lower_ = 2.0;
        
        vector<double> source_distribution_;
        
        vector<Weight_Window> weight_windows_;
        
        Bank bank_;
        Data data_;
        Estimators estimators_;
        Mesh mesh_;
        
        std::mt19937 generator_;
        std::uniform_real_distribution<double> distribution_;
        
        inline double get_rand()
        {
            return distribution_(generator_);
        }

        inline double sample_mfp()
        {
            return - log(get_rand());
        }

        inline double &source_distribution(unsigned cell,
                                          unsigned group)
        {
            return source_distribution_[group + data_.number_of_groups() * cell];
        }

        inline Weight_Window &weight_windows(unsigned cell,
                                             unsigned group)
        {
            return weight_windows_[group + data_.number_of_groups() * cell];
        }
        
        void get_source_particle();
        
        double get_position(unsigned cell);

        void compute_source_distribution();

        void simulate_collision(Particle &particle);

        double get_isotropic_angle();

        double get_anisotropic_angle(unsigned &to_group,
                                     Particle &particle);
        
        int get_cell(double position);

        void stream_left(double &mfp,
                         Particle &particle);

        void stream_right(double &mfp,
                          Particle &particle);

        double distance_to_boundary(unsigned cell,
                                    double position,
                                    double direction);

        void check_weight_windows(Particle &particle);
        
        void initialize_weight_windows();

    public:
    
        Monte_Carlo(unsigned &number_of_cells,
                    unsigned &number_of_groups,
                    unsigned &number_of_scattering_moments,
                    unsigned &number_of_histories,
                    double &side_length,
                    vector<double> &internal_source,
                    vector<double> &boundary_sources,
                    vector<double> &sigma_t,
                    vector<double> &sigma_s,
                    vector<double> &nu,
                    vector<double> &sigma_f,
                    vector<double> &chi,
                    vector<string> &boundary_conditions,
                    bool implicit_capture = true);

        void initialize_weight_windows(vector<double> &phi_adjoint, double max_split = 1.2);
        
        void solve();

        void print_scalar_flux()
        {
            cout << "Monte Carlo";
            cout << "\tTotal particles: " << bank_.total_particles();
            cout << "\tNum histories: " << number_of_histories_;
            cout << "\tImplicit capture: " << implicit_capture_;
            cout << "\tInit with adjoint: " << init_with_adjoint_;
            if (init_with_adjoint_)
            {
                cout << "\tMax split: " << max_splitting_;
            }
            cout << "\tRuntime: " << runtime_;

            cout << endl;
            
            estimators_.print();
        }

    };
}

#endif
