#!/bin/bash --login

#PBS -N domcfg
#PBS -l walltime=00:03:00
#PBS -j oe
#PBS -q normal
#PBS -l select=1
#PBS -P jmmp

#module add cp2k
#module unload cray-snplauncher

  export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)
  export OMP_NUM_THREADS=1
  cd $PBS_O_WORKDIR
#
  #cd /projects/jomp/pmathiot/NEMO_trunk_OVF/tools/DOMAINcfg/eORCA025
  #echo Directory is `pwd`

  OCORES=32
  O_PER_NODE=32

  #OCORES=8
  #O_PER_NODE=8

  #OCORES=1
  #O_PER_NODE=1

  echo time aprun -b  -n $OCORES -N $O_PER_NODE make_domain_cfg.exe
  time aprun -b  -n $OCORES -N $O_PER_NODE ./make_domain_cfg.exe

