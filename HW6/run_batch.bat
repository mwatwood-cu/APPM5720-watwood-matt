#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --time=00:40:00
#SBATCH --partition=shas
#SBATCH --output=sample-%j.out

module purge
module load gcc
module load impi


echo "== This is the scripting step! =="
mpirun -n 16 ./we_mpi.x 
echo "== End of Job =="
