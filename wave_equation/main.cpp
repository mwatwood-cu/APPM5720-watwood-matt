#include <iostream>
#include "WaveSimulation.h"

using namespace std;

int main(int argc, char* argv[])
{
    int grid_points = atoi(argv[1]);
    WaveSimulation w(grid_points, grid_points, 0, 1, 0, 1);
    w.Verbose_Test();

    //w.Calculate_Laplacian();
    //w.print_all();

    return 0;
}
