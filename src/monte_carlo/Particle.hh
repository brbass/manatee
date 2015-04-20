#ifndef Particle_hh
#define Particle_hh

namespace monte_carlo_ns
{
    class Particle
    {
    private:

        bool is_alive_;
        bool in_problem_;
        
        unsigned cell_;
        unsigned group_;
        double weight_;
        double position_;
        double angle_;
        
    public:
        
        Particle(unsigned cell,
                 unsigned group,
                 double weight,
                 double position,
                 double angle);
        
        inline bool &is_alive()
        {
            return is_alive_;
        }
        
        inline bool &in_problem()
        {
            return in_problem_;
        }
        
        inline unsigned &cell()
        {
            return cell_;
        }
        
        inline unsigned &group()
        {
            return group_;
        }
        
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
