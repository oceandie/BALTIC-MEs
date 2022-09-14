#!/bin/bash

#
# setup script for hpg tests on bi
#

# exit when any command fails
set -e

# input handling
if [[ $# -ne 2 ]]; then
    echo "usage: ./create_domaincfg.sh <outputdir> <bathymetry>"
    echo "exiting"
    exit 1
fi

outputdir=$1
bathymetry=`realpath $2`
files_dir=files
ROOTDIR=$PWD

mkdir -p $outputdir

make_domain_cfg=/nobackup/smhid20/users/sm_erimu/NEMO/BALTIC-MEs/src/f90/NEMO_4.0.4_hpge_ovf/tools/DOMAINcfg/make_domain_cfg.exe
ln -sfv $make_domain_cfg $outputdir/make_domain_cfg.exe
ln -sfv $bathymetry $outputdir/bathy_meter.nc

cp -v $files_dir/namelist* $outputdir/.

cd $outputdir

if [[ "$HOSTNAME" =~ n[0-9].* ]]; then
    echo "on node" $BASH_REMATCH
    module purge
    module load buildenv-intel/2018a-eb
    time mpprun -np 1 ./make_domain_cfg.exe
else
    echo "ERROR: not on a node, exiting"
    exit 1
fi

exit 0
