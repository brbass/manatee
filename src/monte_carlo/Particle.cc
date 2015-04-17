#include "Monte_Carlo.hh"

namespace monte_carlo_ns
{
    Particle::
    Particle(unsigned group,
             double weight,
             double position,
             double angle):
        group_(group),
        weight_(weight),
        position_(position),
        angle_(angle)
    {
    }
}
