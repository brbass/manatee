#ifndef Finite_Element_hh
#define Finite_Element_hh

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
        double coeff[2][2] = {{1/2, -1/2},
                              {1/2, 1/2}};
        
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

        inline double stiffness(unsigned kb, unsigned kw, unsigned kdb = 0, unsigned kdw = 0)
        {
            // kb: which basis function
            // kw: which weight function
            // kd: derivative of basis function
            // return integral of the product

            unsigned kd = kdb * kdb + kdw * kdw;
            
            // if (kdb==0 && kdw==0)
            //     kd = 0;
            // if (kdb==0 && kdw==1)
            //     kd = 1;
            // if (kdb==1 && kdw==1)
            //     kd = 2;
            // if (kdb==2 && kdw==0)
            //     kd = 4;
            // if (kdb==2 && kdw==1)
            //     kd = 5;
            // if (kdb==2 && kdw==2)
            //     kd = 8
            
            // linear
            switch (kd)
            {
            case 0:
                return 2 * coeff[kb][0] * coeff[kw][0]
                    + 2/3 * coeff[kb][1] * coeff[kw][1];
            case 1:
                if (kdb > kdw)
                    return 2 * coeff[kb][1] * coeff[kw][0];
                else
                    return 2 * coeff[kb][0] * coeff[kw][1];
            case 2:
                return 2 * coeff[kb][1] * coeff[kw][1];
            default:
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
