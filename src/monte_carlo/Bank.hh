#ifndef Bank_hh
#define Bank_hh

#include <stack>

#include "Particle.hh"

namespace monte_carlo_ns
{
    using std::stack;

    class Bank
    {
    private:
    
        stack<Particle> particle_bank_;
        
        long long total_particles_ = 0;
        
    public:

        Bank();
    
        inline bool is_empty()
        {
            return particle_bank_.empty();
        }

        inline bool particles_remain()
        {
            return !particle_bank_.empty();
        }

        inline unsigned number_of_particles()
        {
            return particle_bank_.size();
        }

        long long total_particles()
        {
            return total_particles_;
        }
        
        void add_particle(Particle &part)
        {
            particle_bank_.push(part);

            total_particles_ += 1;
        }
    
        void add_particle(unsigned cell,
                          unsigned group,
                          double weight,
                          double position,
                          double angle)
        {
            particle_bank_.emplace(cell,
                                   group,
                                   weight,
                                   position,
                                   angle);
            
            total_particles_ += 1;
        }
    
        inline Particle get_particle()
        {
            Particle particle(particle_bank_.top());
            
            particle_bank_.pop();
                
            return particle;
        }
    };
}

#endif

