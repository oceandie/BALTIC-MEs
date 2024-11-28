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
coordinates=`realpath $files_dir/coordinates.nc`

mkdir -p $outputdir

# make_domain_cfg=`realpath ../src/f90/NEMO_4.0.4_hpge_ovf/tools/DOMAINcfg/make_domain_cfg.exe`
make_domain_cfg=`realpath ../src/f90/NEMO_repo/tools/DOMAINcfg/make_domain_cfg.exe`

ln -sfv $make_domain_cfg $outputdir/make_domain_cfg.exe
ln -sfv $bathymetry $outputdir/bathy_meter.nc
ln -sfv $coordinates $outputdir/coordinates.nc

cp -v $files_dir/namelist* $outputdir/.

cd $outputdir

if [[ "$HOSTNAME" =~ n[0-9].* ]]; then
    echo "on a node" $BASH_REMATCH
    time srun ./make_domain_cfg.exe
else
    echo "ERROR not on a node, exiting"
    exit 1
fi
