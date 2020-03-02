#include "BoundaryCondition.h"

BoundaryCondition::BoundaryCondition(){
     boundary_style = 1;
     function_style = 3;
     set_function_bounds(0,0,0);
};


BoundaryCondition::BoundaryCondition(int b_style){
     if(b_style>3 || b_style<1) {
         cout<<"Invalid Boundary Style"<<endl;
         boundary_style = 1;
     }
     else {
         boundary_style = b_style;
     }
     function_style = 3; 
     set_function_bounds(0,0,0);
}

BoundaryCondition::BoundaryCondition(int b_style, int f_style){
     if(b_style>3 || b_style<1) {
         cout<<"Invalid Boundary Style"<<endl;
         boundary_style = 1;
     }
     else {
         boundary_style = b_style;
     }

     if(f_style>3 || f_style<1) {
         cout<<"Invalid Boundary Style"<<endl;
         function_style =3;
         set_function_bounds(0,0,0);   
     }
     else {
         function_style = f_style;
     }
}


double BoundaryCondition::get_boundary_value(double x, double y, double t){
   switch(function_style)
   {
       case 1:
       {
           return A*sin(omega*(t+t+0));
       }
       case 2:
       { 
           double sum = 0;
           for(int i=1; i<=c_values.m_ni; i++){
               sum += c_values(i)*pow(t,i-1);
           }
           return sum;
       }
       case 3:
       {
            return A*exp(-pow(sigma*(t-t_0),2));
       }
       default:
       {
            cout<<"Not a valid function type set "<<function_style<<endl;
            return 0;
       }
   }
}
void BoundaryCondition::set_function_style(int f_style){
    function_style = f_style;
}
void BoundaryCondition::set_boundary_style(int b_style){
    boundary_style = b_style;
}
void BoundaryCondition::set_function_bounds(double A_tmp, double t_0_tmp){
    A = A_tmp;
    t_0 = t_0_tmp;
}
void BoundaryCondition::set_function_bounds(Darray1& c_values_tmp){
    c_values.copy(c_values_tmp);
}
void BoundaryCondition::set_function_bounds(double A_tmp, double sigma_tmp, double t_0_tmp){
    A = A_tmp;
    sigma = sigma_tmp;
    omega = sigma_tmp;
    t_0 = t_0_tmp;
}
