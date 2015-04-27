#include "Estimators.hh"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>

namespace monte_carlo_ns
{
    using namespace std;
    
    Estimators::
    Estimators(unsigned number_of_cells,
               unsigned number_of_surfaces,
               unsigned number_of_groups,
               bool track_variance):
        number_of_cells_(number_of_cells),
        number_of_surfaces_(number_of_surfaces),
        number_of_groups_(number_of_groups),
        ongoing_history_(false),
        number_of_histories_(0)
    {
        current_.size_ = number_of_surfaces_ * number_of_groups_;
        surface_flux_.size_ = number_of_surfaces_ * number_of_groups_;
        volume_flux_.size_ = number_of_cells_ * number_of_groups_;
        collision_flux_.size_ = number_of_cells_ * number_of_groups_;
        
        zero_total_vectors();
    }
    
    void Estimators::
    begin_history()
    {
        if (ongoing_history_)
        {
            cerr << "history already started" << endl;
        }
        else
        {
            ongoing_history_ = true;
            
            zero_history_vectors();
        }
    }
    
    void Estimators::
    end_history()
    {
        if (!ongoing_history_)
        {
            cerr << "history not yet begun" << endl;
        }
        else
        {
            ongoing_history_ = false;
            number_of_histories_ += 1;
            
            sum_history_vectors();
        }
    }
    
    void Estimators::
    zero_history_vectors()
    {
        for (unsigned j = 0; j < estimators_.size(); ++j)
        {
            estimators_[j]->history_.assign(estimators_[j]->size_, 0.0);
        }
    }
    
    void Estimators::
    zero_total_vectors()
    {
        for (unsigned j = 0; j < estimators_.size(); ++j)
        {
            estimators_[j]->total_.assign(estimators_[j]->size_, 0.0);
            estimators_[j]->total_squared_.assign(estimators_[j]->size_, 0.0);
        }
    }
    
    void Estimators::
    sum_history_vectors()
    {
        for (unsigned j = 0; j < estimators_.size(); ++j)
        {
            for (unsigned i = 0; i < estimators_[j]->size_; ++i)
            {
                estimators_[j]->total_[i] += estimators_[j]->history_[i];
                estimators_[j]->total_squared_[i] += pow(estimators_[j]->history_[i], 2);
            }
        }
    }

    void Estimators::
    normalize(double &total_source,
              unsigned &number_of_histories,
              vector<double> &cell_volume)
    {
        for (unsigned i = 0; i < number_of_surfaces_; ++i)
        {
            for (unsigned g = 0; g < number_of_groups_; ++g)
            {
                current_.total_[g + number_of_groups_ * i] *= total_source / number_of_histories;
                current_.total_squared_[g + number_of_groups_ * i] *= total_source / number_of_histories;

                surface_flux_.total_[g + number_of_groups_ * i] *= total_source / number_of_histories;
                surface_flux_.total_squared_[g + number_of_groups_ * i] *= total_source / number_of_histories;
            }
        }

        for (unsigned i = 0; i < number_of_cells_; ++i)
        {
            for (unsigned g = 0; g < number_of_groups_; ++g)
            {
                volume_flux_.total_[g + number_of_groups_ * i] *= total_source / (number_of_histories * cell_volume[i]);
                volume_flux_.total_squared_[g + number_of_groups_ * i] *= total_source / (number_of_histories * cell_volume[i]);
                
                collision_flux_.total_[g + number_of_groups_ * i] *= total_source / (number_of_histories * cell_volume[i]);
                collision_flux_.total_squared_[g + number_of_groups_ * i] *= total_source / (number_of_histories * cell_volume[i]);
            }
        }
    }
    
    void Estimators::
    add_track(unsigned cell,
              unsigned group,
              double weight,
              double length)
    {
        // if (weight * length < 0)
        // {
        //     cout << "track: " << length << endl;
        // }
        
        volume_flux_.history_[group + number_of_groups_ * cell] += weight * length;
    }

    void Estimators::
    add_collision(unsigned cell,
                  unsigned group,
                  double weight,
                  double sigma_t)
    {
        collision_flux_.history_[group + number_of_groups_ * cell] += weight / sigma_t;
    }
    
    void Estimators::
    add_crossing(unsigned surface,
                 unsigned group,
                 double weight,
                 double ordinate)
    {
        current_.history_[group + number_of_groups_ * surface] += weight;
        surface_flux_.total_[group + number_of_groups_ * surface] += weight / abs(ordinate);
    }

    void Estimators::
    print()
    {
        int width = 16;
        
        cout << "Estimator values" << endl << endl;
        
        cout << setw(width) << "cell/edge" << setw(width) << "group";
        for (unsigned j = 0; j < estimator_names_.size(); ++j)
        {
            cout << setw(width) << estimator_names_[j];
        }
        cout << endl;
        
        for (unsigned g = 0; g < number_of_groups_; ++g)
        {
            for (unsigned i = 0; i < number_of_surfaces_; ++i)
            {
                cout << setw(width) << i << setw(width) << g;
                cout << setw(width) << current_.total_[g + number_of_groups_ * i];
                cout << setw(width) << surface_flux_.total_[g + number_of_groups_ * i];
                if (i < number_of_cells_)
                {
                    cout << setw(width) << volume_flux_.total_[g + number_of_groups_ * i];
                    cout << setw(width) << collision_flux_.total_[g + number_of_groups_ * i];
                }

                cout << setw(width) << current_.total_squared_[g + number_of_groups_ * i] - pow(current_.total_[g + number_of_groups_ * i], 2);
                cout << setw(width) << surface_flux_.total_squared_[g + number_of_groups_ * i] - pow(surface_flux_.total_[g + number_of_groups_ * i], 2);
                if (i < number_of_cells_)
                {
                    cout << setw(width) << volume_flux_.total_squared_[g + number_of_groups_ * i] - pow(volume_flux_.total_[g + number_of_groups_ * i], 2);
                    cout << setw(width) << collision_flux_.total_squared_[g + number_of_groups_ * i] - pow(collision_flux_.total_[g + number_of_groups_ * i], 2);
                }
                
                cout << endl;
            }
        }
    }
}
