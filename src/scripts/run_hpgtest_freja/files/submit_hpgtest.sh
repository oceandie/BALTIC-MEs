#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH -N 4
#SBATCH -n 252
#SBATCH -J HPG_TEST
#SBATCH --qos low

module purge
module load buildenv-intel/2023a-eb
module load netCDF-HDF5/4.9.2-1.12.2-hpc1

time srun --multi-prog ./cpu_mapping_freja.conf
