#include "Monte_Carlo.hh"

#include <cmath>
#include <ctime>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include <Teuchos_Time.hpp>

#include "Data.hh"
#include "Mesh.hh"
#include "Particle.hh"

namespace monte_carlo_ns
{
    using namespace std;
    using namespace data_ns;
    using namespace mesh_ns;
    
    Monte_Carlo::
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
                bool implicit_capture):
        data_(number_of_cells,
              number_of_groups,
              number_of_scattering_moments,
              internal_source,
              boundary_sources,
              sigma_t,
              sigma_s,
              nu,
              sigma_f,
              chi,
              boundary_conditions),
        estimators_(number_of_cells,
                    number_of_cells + 1,
                    number_of_groups,
                    true),
        mesh_(number_of_cells,
              side_length),
        number_of_histories_(number_of_histories),
        implicit_capture_(implicit_capture),
	distribution_(0.0, 1.0)
    {
        data_.check();
        mesh_.check();
        
        compute_source_distribution();
        
        initialize_weight_windows();

        generator_.seed(time(nullptr));
        
        // check whether condition sigma_s0 >= 3 * sigma_s1 and sigma_s1 > 0
        // if (data_.number_of_scattering_moments() == 2)
        // {
        //     for (unsigned i = 0; i < mesh_.number_of_cells(); ++i)
        //     {
        //         for (unsigned g = 0; g < data_.number_of_groups(); ++g)
        //         {
        //         }
        //     }
        // }
        
        // add additional checks in Data for this kind of thing, and for things like the sums of cross sections and the total.

        // potential for infinite loops in while statements and recursive function calls should be addressed
    }

    void Monte_Carlo::
    compute_source_distribution()
    {
        source_distribution_.assign(mesh_.number_of_cells() * data_.number_of_groups(), 0);
        
        total_source_ = 0;
        
        for (unsigned i = 0; i < mesh_.number_of_cells(); ++i)
        {
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                total_source_ += mesh_.cell_length(i) * data_.internal_source(i, g);
            }
        }
        
        for (unsigned i = 0; i < mesh_.number_of_cells(); ++i)
        {
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                source_distribution(i, g) = mesh_.cell_length(i) * data_.internal_source(i, g) / total_source_;
            }
        }
    }
    
    void Monte_Carlo::
    get_source_particle()
    {
        double source_sum = 0;
        double source_rand = get_rand();
        
        for (unsigned i = 0; i < mesh_.number_of_cells(); ++i)
        {
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                source_sum += source_distribution(i, g);
                
                if (source_rand < source_sum)
                {
                    double position = get_position(i);
                    
                    bank_.add_particle(i,
                                       g,
                                       birth_weight_,
                                       position,
                                       get_isotropic_angle());
                    
                    //total_birth_weight_ += weight_windows(i, g).ideal;
                    // if (position < mesh_.cell_edge_position(i, 0) || position > mesh_.cell_edge_position(i, 1))
                    // {
                    //     cout << "source  ";
                    //     cout << "  left: " << mesh_.cell_edge_position(i, 0);
                    //     cout << "  right: " << mesh_.cell_edge_position(i, 1);
                    //     cout << "  position: " << position << endl;
                    // }
                    
                    return;
                }
            }
        }
    }
    
    double Monte_Carlo::
    get_position(unsigned cell)
    {
        return mesh_.cell_center_position(cell) + mesh_.cell_length(cell) * (get_rand() - 0.5);
    }
    
    void Monte_Carlo::
    solve()
    {
        Teuchos::Time timer("solve");
        timer.start();
        
        for (unsigned h = 0; h < number_of_histories_; ++h)
        {
            get_source_particle();

            estimators_.begin_history();

            while(bank_.particles_remain())
            {
                
                // if (bank_.number_of_particles() % int(1e4) == 0 && bank_.number_of_particles() > 1e4)
                // {
                //     cout << "number of particles in bank: ";
                //     cout << bank_.number_of_particles() << endl;
                // }
                
                Particle particle = bank_.get_particle();

                // if (particle.position() < mesh_.cell_edge_position(particle.cell(), 0) || particle.position() > mesh_.cell_edge_position(particle.cell(), 1))
                // {
                //     cout << "angle: " << particle.angle();
                //     cout << "left: " << mesh_.cell_edge_position(particle.cell(), 0);
                //     cout << "  right: " << mesh_.cell_edge_position(particle.cell(), 1);
                //     cout << "  position: " << particle.position() << endl;
                // }
                
                check_weight_windows(particle); 

                while (particle.is_alive())
                {
                    double mfp = sample_mfp();
                    
                    while (mfp > 0) // distance remains in sample
                    {
                        if (particle.angle() > 0)
                        {
                            stream_right(mfp,
                                         particle);
                        }
                        else
                        {
                            stream_left(mfp,
                                        particle);
                        }

                        if (particle.is_alive())
                        {
                            check_weight_windows(particle);
                        }
                        else
                        {
                            break;
                        }
                    } 
                    
                    if (particle.is_alive())
                    {
                        simulate_collision(particle);
                        
                        check_weight_windows(particle);
                    }
                }
            }
            estimators_.end_history();
        }
        
        vector<double> cell_volume(mesh_.number_of_cells());
        
        for (unsigned i = 0; i < mesh_.number_of_cells(); ++i)
        {
            cell_volume[i] = mesh_.cell_length(i);
        }
        
        estimators_.normalize(total_source_,
                              number_of_histories_,
                              cell_volume);
        
        runtime_ = timer.stop();
    }
    
    void Monte_Carlo::
    simulate_collision(Particle &particle)
    {
        estimators_.add_collision(particle.cell(),
                                  particle.group(),
                                  particle.weight(),
                                  data_.sigma_t(particle.cell(), particle.group()));
        
        double cross_rand = 0;
        
        if (implicit_capture_)
        {
            particle.weight() *= 1 - data_.sigma_a(particle.cell(), particle.group()) / data_.sigma_t(particle.cell(), particle.group());
            
            cross_rand = get_rand() * (data_.sigma_t(particle.cell(), particle.group()) - data_.sigma_a(particle.cell(), particle.group()));
        }
        else
        {
            cross_rand = get_rand() * data_.sigma_t(particle.cell(), particle.group());
        }
        
        double cross_sum = 0;
        
        // scattering
        for (unsigned g = 0; g < data_.number_of_groups(); ++g)
        {
            cross_sum += data_.sigma_s(particle.cell(), particle.group(), g);
            
            if (cross_rand < cross_sum)
            {
                if (data_.number_of_scattering_moments() == 1)
                {
                    particle.angle() = get_isotropic_angle();
                }
                else if (data_.number_of_scattering_moments() == 2)
                {
                    particle.angle() = get_anisotropic_angle(g,
                                                             particle);
                }
                else
                {
                    cerr << "only first moment of scattering supported" << endl;
                }
                
                particle.group() = g;
                
                return;
            }
        }
        
        // fission
        cross_sum += data_.sigma_f(particle.cell(), particle.group());
        
        if (cross_rand < cross_sum)
        {
            particle.is_alive() = false;
            
            unsigned num_fission_particles = floor(data_.nu(particle.cell(), particle.group()) + get_rand());
            
            for (unsigned p = 0; p < num_fission_particles; ++p)
            {
                double chi_rand = get_rand();
                double chi_sum = 0;
                
                for (unsigned g = 0; g < data_.number_of_groups(); ++g)
                {
                    chi_sum += data_.chi(particle.cell(), g);
                    
                    if (chi_rand < chi_sum)
                    {
                        bank_.add_particle(particle.cell(),
                                           g,
                                           particle.weight(),
                                           particle.position(),
                                           get_isotropic_angle());
                        
                        break;
                    }
                }
            }
        }
        
        // absorption
        cross_sum += data_.sigma_a(particle.cell(), particle.group());
        
        if (cross_rand < cross_sum)
        {
            particle.is_alive() = false;
            
            return;
        }

        // cross sections not equal to total
        cerr << "cross sections do not sum to total" << endl;
    }

    double Monte_Carlo::
    get_isotropic_angle()
    {
        return 2 * get_rand() - 1;
    }

    double Monte_Carlo::
    get_anisotropic_angle(unsigned &to_group,
                          Particle &particle)
    {
        if (data_.sigma_s(particle.cell(), particle.group(), to_group, 1) > 0)
        {
            double s1 = pow(data_.sigma_s(particle.cell(), particle.group(), to_group, 0) - 3 * particle.angle() * data_.sigma_s(particle.cell(), particle.group(), to_group, 1), 2);
            double s2 = 12 * particle.angle() * data_.sigma_s(particle.cell(), particle.group(), to_group, 0) * data_.sigma_s(particle.cell(), particle.group(), to_group, 1) * get_rand();
            double s3 = 3 * particle.angle() * data_.sigma_s(particle.cell(), particle.group(), to_group, 1);
            
            return (- data_.sigma_s(particle.cell(), particle.group(), to_group, 0) + sqrt(s1 + s2)) / s3;
        }
        else
        {
            return get_isotropic_angle();
        }
    }
    
    int Monte_Carlo::
    get_cell(double position)
    {
        for (unsigned i = 0; i < mesh_.number_of_cells(); ++i)
        {
            if (position < mesh_.cell_edge_position(i, 1))
            {
                return i;
            }
        }

        return -1;
    }

    void Monte_Carlo::
    stream_left(double &mfp,
                Particle &particle)
    {
        double distance = mfp / data_.sigma_t(particle.cell(), particle.group());
        double boundary_distance = distance_to_boundary(particle.cell(), particle.position(), particle.angle());

        if (distance > boundary_distance) // cross cell boundary
        {
            mfp -= boundary_distance * data_.sigma_t(particle.cell(), particle.group());
            
            estimators_.add_track(particle.cell(),
                                  particle.group(),
                                  particle.weight(),
                                  boundary_distance);
            
            estimators_.add_crossing(particle.cell(),
                                     particle.group(),
                                     particle.weight(),
                                     particle.angle());
            
            particle.position() = mesh_.cell_edge_position(particle.cell(), 0);
            
            if (particle.cell() == 0) // at left boundary
            {
                if (data_.boundary_condition(0).compare("reflected") == 0)
                {
                    particle.angle() = - particle.angle();
                    
                    estimators_.add_crossing(particle.cell(),
                                             particle.group(),
                                             particle.weight(),
                                             particle.angle());
                }
                else if (data_.boundary_condition(0).compare("vacuum") == 0)
                {
                    mfp = 0;
                    
                    particle.is_alive() = false;
                }
                else
                {
                    cerr << "boundary condition does not exist" << endl;
                }
            }
            else // at another cell boundary
            {
                particle.cell() -= 1;
            }
        }
        else
        {
            estimators_.add_track(particle.cell(),
                                  particle.group(),
                                  particle.weight(),
                                  distance);

            mfp = 0;
            
            particle.position() += distance * particle.angle();
        }
    }

    void Monte_Carlo::
    stream_right(double &mfp,
                 Particle &particle)
    {
        double distance = mfp / data_.sigma_t(particle.cell(), particle.group());
        double boundary_distance = distance_to_boundary(particle.cell(), particle.position(), particle.angle());
        
        if (distance > boundary_distance) // cross cell boundary
        {
            mfp -= boundary_distance * data_.sigma_t(particle.cell(), particle.group());
            
            estimators_.add_track(particle.cell(),
                                  particle.group(),
                                  particle.weight(),
                                  boundary_distance);
            
            estimators_.add_crossing(particle.cell() + 1,
                                     particle.group(),
                                     particle.weight(),
                                     particle.angle());
            
            particle.position() = mesh_.cell_edge_position(particle.cell(), 1);
            
            if (particle.cell() == mesh_.number_of_cells() - 1) // at left boundary
            {
                if (data_.boundary_condition(1).compare("reflected") == 0)
                {
                    particle.angle() = - particle.angle();
                    
                    estimators_.add_crossing(particle.cell() + 1,
                                             particle.group(),
                                             particle.weight(),
                                             particle.angle());
                }
                else if (data_.boundary_condition(1).compare("vacuum") == 0)
                {
                    mfp = 0;
                    
                    particle.is_alive() = false;
                }
                else
                {
                    cerr << "boundary condition does not exist" << endl;
                }
            }
            else // at another cell boundary
            {
                particle.cell() += 1;
            }
        }
        else
        {
            estimators_.add_track(particle.cell(),
                                  particle.group(),
                                  particle.weight(),
                                  distance);
            
            mfp = 0;
            
            particle.position() += distance * particle.angle();
        }
    }

    double Monte_Carlo::
    distance_to_boundary(unsigned cell,
                         double position,
                         double angle)
    {
        if (angle > 0)
        {
            return (mesh_.cell_edge_position(cell, 1) - position) / angle;
        }
        else
        {
            return (position - mesh_.cell_edge_position(cell, 0)) / abs(angle);
        }
    }
    
    void Monte_Carlo::
    check_weight_windows(Particle &particle)
    {
        if (particle.weight() < weight_windows(particle.cell(), particle.group()).lower)
        {
            double survival_probability = particle.weight() / weight_windows(particle.cell(), particle.group()).ideal;
            
            if (get_rand() < survival_probability)
            {
                particle.weight() = weight_windows(particle.cell(), particle.group()).ideal;
            }
            else
            {
                particle.is_alive() = false;
            }
        }
        else if (particle.weight() > weight_windows(particle.cell(), particle.group()).upper) // assumes ideal <= 2 * upper
        {
            unsigned particles_to_create = floor(particle.weight() / weight_windows(particle.cell(), particle.group()).ideal) - 1;
            
            particle.weight() /= (particles_to_create + 1);
            
            for (unsigned p = 0; p < particles_to_create; ++p)
            {
                bank_.add_particle(particle.cell(),
                                   particle.group(),
                                   particle.weight(),
                                   particle.position(),
                                   particle.angle());
            }
        }
    }

    void Monte_Carlo::
    initialize_weight_windows()
    {
        weight_windows_.resize(mesh_.number_of_cells() * data_.number_of_groups());
    }
    
    void Monte_Carlo::
    initialize_weight_windows(vector<double> &phi_adjoint)
    {
        init_with_adjoint_ = true;
        
        weight_windows_.resize(mesh_.number_of_cells() * data_.number_of_groups());
        
        if (phi_adjoint.size() != mesh_.number_of_cells() * mesh_.number_of_nodes() * data_.number_of_groups())
        {
            cerr << "phi adjoint is incorrect size" << endl;
            return;
        }
        vector<double> phi(mesh_.number_of_cells() * data_.number_of_groups());

        // calculate average phi for a cell
        for (unsigned i = 0; i < mesh_.number_of_cells(); ++i)
        {
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                double sum = 0;
                    
                for (unsigned n = 0; n < mesh_.number_of_nodes(); ++n)
                {
                    unsigned k = n + mesh_.number_of_nodes() * (g + data_.number_of_groups() * i);
                        
                    sum += abs(phi_adjoint[k]);
                }
                    
                unsigned k = g + data_.number_of_groups() * i;
                    
                phi[k] = sum / mesh_.number_of_nodes();
            }
        }

        // calculate the average value of phi
        double response = 0;
        double phi_average = 0;
        for (unsigned i = 0; i < mesh_.number_of_cells(); ++i)
        {
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                unsigned k = g + data_.number_of_groups() * i;

                // response += phi[k] * mesh_.cell_length(i) * data_.sigma_t(i, g);
                phi_average += phi[k] * mesh_.cell_length(i);
            }
        }
        phi_average /= mesh_.side_length();

        // use these values to calculate the weight windows
        for (unsigned i = 0; i < mesh_.number_of_cells(); ++i)
        {
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                unsigned k = g + data_.number_of_groups() * i;
                
                weight_windows(i, g).lower = phi_average / phi[k];

                weight_windows(i, g).upper = weight_windows(i, g).lower * upper_over_lower_;
                
                weight_windows(i, g).ideal = weight_windows(i, g).lower * ideal_over_lower_;
            }
        }

        // calculate the average weight window for each group at the points where particles are born
        double average_source_weight = 0;
        double average_source = 0;
        
        for (unsigned i = 0; i < mesh_.number_of_cells(); ++i)
        {
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                
                average_source_weight += weight_windows(i, g).ideal * mesh_.cell_length(i) * data_.internal_source(i, g);
                average_source += mesh_.cell_length(i) * data_.internal_source(i, g);
            }
        }
        
        average_source_weight /= average_source;
        
        for (unsigned i = 0; i < mesh_.number_of_cells(); ++i)
        {
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                weight_windows(i, g).ideal /= average_source_weight;
                weight_windows(i, g).lower = weight_windows(i, g).ideal / ideal_over_lower_;
                weight_windows(i, g).upper = weight_windows(i, g).lower * upper_over_lower_;
            }
        }
    }
}
