#include "Ordinates.hh"

#include "Gauss_Legendre.hh"

namespace transport_ns
{
    Ordinates::
    Ordinates(unsigned number_of_ordinates):
        number_of_ordinates_(number_of_ordinates)
    {
        gauss_legendre_vec(number_of_ordinates_, ordinates_, weights_);

        alpha_.resize(number_of_ordinates_, 0);
        alpha_half_.resize(number_of_ordinates_, 0);
        
        alpha_half_[0] = -ordinates_[0] * weights_[0];
        alpha_[0] = alpha_half_[0];
        
        for (unsigned o = 1; o < number_of_ordinates_; ++o)
        {
            alpha_half_[o] = alpha_half_[o-1] - ordinates_[o] * weights_[o];
            alpha_[o] = alpha_half_[o] + alpha_half_[o-1];
        }
    }
}
