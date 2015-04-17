#ifndef Transport_Model_hh
#define Transport_Model_hh

#include <iostream>

namespace transport_ns
{
    class Transport_Model
    {
        
    protected:

    public:
        
        Transport_Model();
        
        virtual void solve()
        {
            std::cerr << "solve not implemented" << std::endl;
        }
        
        virtual void solve_eigenvalue()
        {
            std::cerr << "solve_eigenvalue not implemented" << std::endl;
        }
        
        virtual void print_scalar_flux()
        {
            std::cerr << "print_scalar_flux not implemented" << std::endl;
        }

        void print_eigenvalue()
        {
            std::cerr << "print_eigenvalue not implemented" << std::endl;
        }

    };
}

#endif
