#ifndef Bank_hh
#define Bank_hh

#include <stack>

#include "Particle.hh"

namespace monte_carlo_ns
{
    using std::stack;
    
private:
    
    stack<Particle> particle_bank_;

    long long num_particles_ = 0;
    
public:

    Bank();
    
    inline bool is_empty()
    {
        return particle_bank_.empty();
    }
    
    void add_particle(Particle &part)
    {
        particle_bank_.push(part);

        num_particles_ += 1;
    }
    
    void add_particle(unsigned group,
                      double weight,
                      double position,
                      double angle)
    {
        particle_bank_.emplace(group,
                               weight,
                               position,
                               angle);

        num_particles_ += 1;
    }
    
    inline Particle get_particle()
    {
        Particle particle(particle_bank_.top());
        particle_bank_.pop();
        
        return particle;
    }
    
}

#endif

