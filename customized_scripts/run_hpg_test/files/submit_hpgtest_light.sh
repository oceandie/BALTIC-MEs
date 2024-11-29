#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH --nodes=4
#SBATCH --ntasks=249
#SBATCH --cpus-per-task=1
#SBATCH -J HPG_TEST

source ~/.bashrc

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

module load buildenv-intel/2023a-eb NCO/5.1.3-hpc1-gcc-2022a-eb \
       netCDF-HDF5/4.9.2-1.12.2-hpc1

printenv > env_vars
module list > modules

#time mpirun --bind-to core -np 64 ./nemo.exe
#time mpirun -np $SLURM_NTASKS ./nemo.exe
time srun --cpu-bind=cores ./nemo.exe

# bash rebuild.sh
