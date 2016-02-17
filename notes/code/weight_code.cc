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
    else if (particle.weight() > weight_windows(particle.cell(), particle.group()).upper)
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
