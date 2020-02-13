#include <iostream>
#include "WaveSimulation.h"

using namespace std;

int main()
{
    WaveSimulation w(1000, 1000, 0, 1, 0, 1);
    w.Verbose_Test();

    //w.Calculate_Laplacian();
    //w.print_all();

    return 0;
}
