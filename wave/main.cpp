#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono> 
#include "silo.h"
#include "Darrays.h"
#include "WaveSimulation.h"
#include "omp.h"

using namespace std::chrono;
using namespace std;

int main()
{
auto start = high_resolution_clock::now();
auto start_loop = high_resolution_clock::now();
auto stop = high_resolution_clock::now();
auto duration = duration_cast<microseconds>(stop - start);
double t;
int initial_data, material_model, forcing, time_source, it;
double energy, energyold;
//Create the sim
WaveSimulation w_sim;
#pragma omp parallel
{
#pragma omp single
{
    start = high_resolution_clock::now();
    t = 0.0;
    
        // Setup of problem
        // These can also be read using cin at run time
    initial_data = 0;
    material_model = 2;
    forcing = 1;
    time_source = 1;
    
    //Setup Boundary Condition style and function.
    w_sim.m_bc.set_boundary_style(1);
    w_sim.m_bc.set_function_style(1);
    w_sim.m_bc.set_function_bounds(0,1,0);

    w_sim.m_check_energy = false;        
    w_sim.m_tend = 10;
    w_sim.set_n_print(40);        
}
        // Call methods that set up problem
    w_sim.set_up_initial_data(initial_data);
    w_sim.set_up_material(material_model);
    w_sim.set_time_source(time_source);
    w_sim.set_up_forcing(forcing);
#pragma omp single
{
    // The timestep of the method must not break the physics
    w_sim.set_up_time_step();
}
   // Make sure boundary conditions are set.
    w_sim.set_boundary_conditions(0);
#pragma omp single
{
    // For checking energy
    energy, energyold = 1.0;
    w_sim.print_solution(0,t);

    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    cout<<"Initialize Took: "<<duration.count()/1000<<" ms"<<endl;
    //Start Timing of Looping
    start_loop = high_resolution_clock::now();
    start = high_resolution_clock::now();
}
} //End Initial parallel region

    for(it = 1; it<= w_sim.m_Nt; it++){
#pragma omp parallel
{
        t = (it-1)*w_sim.m_dt;
        w_sim.set_time(t);
        w_sim.compute_laplace();
        w_sim.advance();
        t = it*w_sim.m_dt;
#pragma omp single
{
        // Print every so often, always print last solution
        if(w_sim.m_check_energy && it%500 == 0){
            stop = high_resolution_clock::now(); 
            duration = duration_cast<microseconds>(stop - start);
            cout<<"1 Loop Took an Average: "<<(double)duration.count()/1000/500<<" ms"<<endl;
            //Starting again
            start = high_resolution_clock::now();
            energy = w_sim.compute_energy();
            cout <<"Energy at time " << setprecision(3) << t << " is: "
                 << setprecision(16) << scientific << right << setw(22) << energy
                 << " diff is: " << setprecision(2)
                 << scientific << right << setw(8)
                 << (energy - energyold)/energy 
                 << "\n";
            energyold = energy;
        }
}
        // Interchange the time levels
        w_sim.swap_times();
            // Set Boundaryy conditions
        w_sim.set_boundary_conditions(t);
}
    }
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop-start_loop);
    cout<<"All Loops Took an Average: "<<(double)duration.count()/1000<<" ms"<<endl; 
    return 0;
}
