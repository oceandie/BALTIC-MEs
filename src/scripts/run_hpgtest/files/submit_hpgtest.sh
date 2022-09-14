#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH -N 37
#SBATCH -n 584
#SBATCH -J HPG_TEST

module purge
module load netCDF/4.4.1.1-HDF5-1.8.19-nsc1-intel-2018a-eb
module load HDF5/1.8.19-nsc1-intel-2018a-eb

np=$SLURM_NTASKS
xios_np=32
#nemo_np=$((np - xios_np))
nemo_np=552

time mpprun -np $nemo_np ./nemo.exe : -np $xios_np ./xios_server.exe

exit
