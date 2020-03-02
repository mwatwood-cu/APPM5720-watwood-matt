#include <iostream>
#include <cmath>
#include "Darrays.h"
#include "BoundaryCondition.h"

using namespace std;

int main(){

    BoundaryCondition aBc(1);
    aBc.set_function_style(2);
    Darray1 cs(1, 2);
    cs(1) = 1;
    cs(2) = 2;
    aBc.set_function_bounds(cs);
    cout<<"Value "<<aBc.get_boundary_value(0,0,1)<<endl;
    return 0;

}

