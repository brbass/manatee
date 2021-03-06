#ifndef Ordinates_hh
#define Ordinates_hh

#include <vector>

#include "Gauss_Legendre.hh"

namespace transport_ns
{
    using std::vector;
    
    class Ordinates
    {
    private:

        unsigned number_of_ordinates_;

        vector<double> ordinates_;
        vector<double> weights_;
        vector<double> alpha_;
        vector<double> alpha_half_;
        
    public:
        
        Ordinates(unsigned number_of_ordinates);

        inline unsigned number_of_ordinates()
        {
            return number_of_ordinates_;
        }

        inline double ordinates(unsigned o)
        {
            return ordinates_[o];
        }

        inline double weights(unsigned o)
        {
            return weights_[o];
        }

        inline double alpha(unsigned o)
        {
            return alpha_[o];
        }

        inline double alpha_half(unsigned o)
        {
            return alpha_half_[o];
        }
    };
}

#endif
