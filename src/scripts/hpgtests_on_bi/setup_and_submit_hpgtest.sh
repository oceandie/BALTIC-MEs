#!/bin/bash

#
# setup script for hpg tests on bi
#

# exit when any command fails
set -e

# input handling
if [[ $# -ne 2 ]]; then
    echo "usage: ./setup_hpgtest.sh <testname> <domcfg>"
    echo "exiting"
    exit
fi

testdir=$1
domcfg=`realpath $2`
template_dir=test_template
ROOTDIR=$PWD

mkdir -p $testdir

nemo=/nobackup/smhid20/users/sm_erimu/NEMO/BALTIC-MEs/src/f90/NEMO_4.0.4_hpge_ovf/cfgs/MESTEST/BLD/bin/nemo.exe
xios=/nobackup/smhid20/users/sm_jongr/NEMO/MODELS/xios-2.5_trunk/bin/xios_server.exe

ln -sfv $nemo $testdir/nemo.exe
ln -sfv $xios $testdir/xios_server.exe

cp -rv $template_dir/submit_hpgtest.sh $testdir/.
cp -v $template_dir/*.xml $testdir/.
cp -v $template_dir/namelist_* $testdir/.

ln -sfv $domcfg $testdir/domain_cfg.nc

cd $testdir
sbatch submit_hpgtest.sh

cd $ROOTDIR
