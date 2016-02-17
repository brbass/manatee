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
