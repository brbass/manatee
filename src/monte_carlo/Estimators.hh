#ifndef Estimators_hh
#define Estimators_hh

#include <iostream>
#include <string>
#include <vector>

namespace monte_carlo_ns
{
    using std::string;
    using std::vector;

    class Estimators
    {
    private:

        bool ongoing_history_;
        unsigned number_of_histories_;
        
        unsigned number_of_cells_;
        unsigned number_of_surfaces_;
        unsigned number_of_groups_;
        
        unsigned number_of_estimators_;
        
        struct Estimator
        {
            unsigned size_;
            vector<double> history_;
            vector<double> total_;
            vector<double> total_squared_;
        };
        
        Estimator current_;
        Estimator surface_flux_;
        Estimator volume_flux_;
        Estimator collision_flux_;
        
        vector<Estimator*> estimators_ = {&current_,
                                          &surface_flux_,
                                          &volume_flux_,
                                          &collision_flux_};
        
        vector<string> estimator_names_ = {"current", "surface_flux", "volume_flux", "collision_flux",
                                           "curr_var", "sur_var", "vol_var", "coll_var",
                                           "curr_err", "sur_err", "vol_err", "coll_err",
                                           "curr_rel", "sur_rel", "vol_rel", "coll_rel",
                                           "curr_fom", "sur_fom", "vol_fom", "coll_fom"};

        void zero_history_vectors();
        void zero_total_vectors();
        void sum_history_vectors();
        
    public:
        
        Estimators(unsigned number_of_cells,
                   unsigned number_of_surfaces,
                   unsigned number_of_groups,
                   bool track_variance);
        
        void begin_history();
        void end_history();
        
        void normalize(double &total_source,
                       unsigned &number_of_histories,
                       vector<double> &cell_volume);
        
        void add_track(unsigned cell,
                       unsigned group,
                       double weight,
                       double length);
        
        void add_collision(unsigned cell,
                           unsigned group,
                           double weight,
                           double sigma_t);
        
        void add_crossing(unsigned surface,
                          unsigned group,
                          double weight,
                          double ordinate);
        
        void print(double runtime);
    };
}

#endif
