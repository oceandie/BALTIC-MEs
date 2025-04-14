#!/bin/bash

#
# setup and run hpg test on Freja
#

# exit when any command fails
set -e

# input handling
if [[ $# -ne 2 ]]; then
    echo "usage: ./run_hpgtest.sh <testname> <domcfg>"
    echo "exiting"
    exit
fi

testdir=$1
domcfg=`realpath $2`
files_dir=files
ROOTDIR=$PWD

mkdir -p $testdir

MODEL=NORDIC_HPGTEST_FREJA_422_fixes
BLDFOLDER=/nobackup/smhid20/users/sm_erimu/NEMO/nemo4/cfgs/${MODEL}/BLD/
nemo=${BLDFOLDER}/bin/nemo.exe
xios=/home/sm_erimu/Projects/XIOS/trunk/seq_build/bin/xios_server.exe

ln -sfv $nemo $testdir/nemo.exe
ln -sfv $xios $testdir/xios_server.exe

cp -v $files_dir/submit_hpgtest.sh $testdir/.
cp -v $files_dir/cpu_mapping_freja.conf $testdir/.
cp -v $files_dir/submit_hpgtest_light.sh $testdir/.
cp -v $files_dir/cpu_mapping_freja_light.conf $testdir/.
cp -v $files_dir/*.xml $testdir/.
cp -v $files_dir/namelist_* $testdir/.
cp -v $files_dir/rebuild.sh $testdir/.

ln -sfv $domcfg $testdir/domain_cfg.nc

cd $testdir
sbatch submit_hpgtest.sh
# sbatch submit_hpgtest_light.sh
