#include  <iostream>
#include <cmath>
#include "Darrays.h"

using namespace std;

class BoundaryCondition
{
    public:
        BoundaryCondition();
        BoundaryCondition(int b_style);
        BoundaryCondition(int b_style, int f_style);
        ~BoundaryCondition(){}
        double get_boundary_value(double x, double y, double t);
        void set_function_style(int f_style);
        void set_function_bounds(double A_tmp, double t_0_tmp);
        void set_function_bounds(Darray1& c_values_tmp);
        void set_function_bounds(double A_tmp, double sigma_tmp, double t_0_tmp);
    private:
        int boundary_style;
        int function_style;
        double A;
        double omega;
        double t_0;
        double sigma;
        Darray1 c_values;
        
};
