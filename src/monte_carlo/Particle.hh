#ifndef Particle_hh
#define Particle_hh

namespace monte_carlo_ns
{
    class Particle
    {
    private:

        unsigned group_;
        double weight_;
        double position_;
        double angle_;
        
    public:
        
        Particle();

        inline double &weight()
        {
            return weight_;
        }
        
        inline double &position()
        {
            return position_;
        }
        
        inline double &angle()
        {
            return angle_;
        }

    };
}

#endif
