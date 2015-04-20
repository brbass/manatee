#include "Particle.hh"

#include "Monte_Carlo.hh"

namespace monte_carlo_ns
{
    Particle::
    Particle(unsigned cell,
             unsigned group,
             double weight,
             double position,
             double angle):
        cell_(cell),
        group_(group),
        weight_(weight),
        position_(position),
        angle_(angle),
        is_alive_(true),
        in_problem_(true)
    {
    }
}
