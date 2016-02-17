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
