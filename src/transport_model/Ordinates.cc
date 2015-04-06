#include "Ordinates.hh"

#include "Gauss_Legendre.hh"

namespace transport_ns
{
    Ordinates::
    Ordinates(unsigned number_of_ordinates):
        number_of_ordinates_(number_of_ordinates)
    {
        gauss_legendre_vec(number_of_ordinates_, ordinates_, weights_);
    }
}
