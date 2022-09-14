#!/bin/bash

#
# setup and run hpg test on bi
#

# exit when any command fails
set -e

# input handling
if [[ $# -ne 2 ]]; then
    echo "usage: ./run_hpgtest.sh <testname> <domcfg>"
    echo "exiting"
    exit 1
fi

testdir=$1
domcfg=`realpath $2`
files_dir=files
ROOTDIR=$PWD

mkdir -p $testdir

nemo=/nobackup/smhid20/users/sm_erimu/NEMO/BALTIC-MEs/src/f90/NEMO_4.0.4_hpge_ovf/cfgs/HPGTEST/BLD/bin/nemo.exe
xios=/nobackup/smhid20/users/sm_jongr/NEMO/MODELS/xios-2.5_trunk/bin/xios_server.exe

ln -sfv $nemo $testdir/nemo.exe
ln -sfv $xios $testdir/xios_server.exe

cp -v $files_dir/submit_hpgtest.sh $testdir/.
cp -v $files_dir/*.xml $testdir/.
cp -v $files_dir/namelist_* $testdir/.

ln -sfv $domcfg $testdir/domain_cfg.nc

cd $testdir
sbatch submit_hpgtest.sh
