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
            Particle particle = bank_.get_particle();
            
            check_weight_windows(particle); 
            
            while (particle.is_alive())
            {
                double mfp = sample_mfp();
                
                while (mfp > 0)
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
        
    estimators_.normalize(total_source_,
                          number_of_histories_,
                          cell_volume_);
    
    runtime_ = timer.stop();
}