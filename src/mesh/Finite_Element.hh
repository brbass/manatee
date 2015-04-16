#ifndef Finite_Element_hh
#define Finite_Element_hh

#include <iostream>

namespace mesh_ns
{
    class Finite_Element
    {
    private:
        
        double cell_center_position_;
        double cell_length_;
        // quadratic
        // double coeff[3][3] = {{0, -1/2, 1/2},
        //                       {1, 0 , -1},
        //                       {0, 1/2, 1/2}};
        double coeff[2][2] = {{1.0 / 2.0, -1.0 / 2.0},
                              {1.0 / 2.0, 1.0 / 2.0}};
        
    public:

        Finite_Element(){}
        
        Finite_Element(double cell_center_position,
                       double cell_length);
        
        inline double cell_center_position()
        {
            return cell_center_position_;
        }
        
        inline double cell_length()
        {
            return cell_length_;
        }

        inline double cell_edge_position(unsigned edge)
        {
            switch(edge)
            {
            case 0:
                return cell_center_position_ - cell_length_ / 2;
            case 1:
                return cell_center_position_ + cell_length_ / 2;
            default:
                std::cerr << "edge not available" << std::endl;
                return 0;
            }
        }
        
        inline double stiffness(unsigned kb, unsigned kw, unsigned kdb = 0, unsigned kdw = 0)
        {
            // kb: which basis function
            // kw: which weight function
            // kd: derivative of basis function
            // return integral of the product
            
            unsigned kd = kdb * kdb + kdw * kdw;
            
            // linear
            switch (kd)
            {
            case 0:
                return 2.0 * coeff[kb][0] * coeff[kw][0]
                    + 2.0/3.0 * coeff[kb][1] * coeff[kw][1];
                
            case 1:
                if (kdb > kdw)
                    return 2.0 * coeff[kb][1] * coeff[kw][0];
                else
                    return 2.0 * coeff[kb][0] * coeff[kw][1];
                
            case 2:
                return 2.0 * coeff[kb][1] * coeff[kw][1];
                
            default:
                std::cout << "stiffness case " << kd << " not available" << std::endl;
                return 0;
            }
            
            // quadratic
            // switch (kd)
            // {
            // case 0:
            //     return 2 * coeff[kb][0] * coeff[kw][0]
            //         + 2/3 * (coeff[kb][1] * coeff[kw][1] + coeff[kb][2] * coeff[kw][0] + coeff[kb][0] * coeff[kw][2])
            //         + 2/5 * coeff[kb][2] * coeff[kw][2];
                
            // case 2:
            //     return 4 * coeff[kb][0] * coeff[kw][2]
            //         + 4/3 * coeff[kb][2] * coeff[kw][2];
            // default:
            //     return 0;
            // }
        }   
    };
}

#endif
