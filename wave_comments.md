# Comments on the wave.f90 file.

Overall, this code follows the same structure given for scientific code. Preamble, Initialize, Main Loop. <br><br>

In the realm of parallelization, we can run this either on a single node or multiple nodes. On multiple nodes, we could break up this domain and give separate nodes separate parts of the domain and connect them using ghosted nodes to handle the communication and updating values along the boundary of each domain section. This would require some coordination from a master node to coordinate I/O and domain decomposition. This would create a tree type of organization where a master node sets up the problem, then hands out tasks, then gathers the results at the end. If the domain is small enough that the communication between nodes would dominate the run time then this can be run on a single node with multiple threads involved. As this is a simpler case I will primarily discuss how to implement that and where speed ups could be added for this.<br><br>

I don't remember the most efficient double looping structure for Fortran (row major vs column major) but for any matrix this should be tested to ensure the optimal memory traversal for the language.To get really saucy, these looping procedures could additionally be broken up to implement vectorization if the hardware being run on can support it. This would involve more particular memory management especiallly memory localization to ensure optimal performance nubmer of accesses.<br><br>

The preamble section creates the needed variables for the simulation. This section seems to be well organized, and comprehensive, but some of the variable names are not particularly clear. For example:<br>
```real(dp) :: tend  =  1.0_dp``` <br>
```real(dp), allocatable, dimension(:,:) :: um,u,up,lap,c2 ```<br>

are ambiguous to me and even nx, ny could be ambiguous to people new to pde solving. There could also be more structural guidance in the code and comments to differentiate the different sections of the code. <br><br>

The Initialize section from lines 27 - 58 takes the created variables and assigns them values. This is setting up our domain, function values within the domain, boundary conditions, and the matrix used as the discretization. As to speed ups this has room for parallelization where each matrix and array calculation like the<br>
``` call set_bc(u,x,y,t,dx,dy,nx,ny) ```<br>
``` call compute_laplace(lap,u,c2,dx,dy,nx,ny)```<br>
calls, can be threaded out using OpenMP or a similar shared memory parallel implementation. This could be done in the domain decomposition example as well as long as the shared memory node had a shared memory form of parallelization on it.<br><br>

In the main loop, we take our initial condition and update it everywhere in the domain, then update our function values after the whole domain is complete, then we move on to the next time step and when the requisite number of steps is complete then we finish that stepping and print out rhe resulting function value.<br>
The primary optimization is to locally (shared memory threads) split the calculations up using OpenMP or something similar and optimizing the compute laplacian. This optimization may leave the laplacian calculation to be done separatly, but could also integrate that into the individual domain steps.<br><br>

Finally in the output, if multiple nodes are used then a parallel I/O like PNetCDF or something similar could be used to minimize the number of node to node communications and dramatically speed up output time.