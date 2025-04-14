#!/bin/bash

module load parallel/20220722-hpc1-gcc-2022a-eb

expname=hpg_test_1d
files=$(ls $expname*_0000.nc)

# create list of basenames
basenames=()
for fl in ${files[@]}; do
    basename=$(echo $fl | sed 's/_0000.nc//')
    basenames+=($basename)
done
initfile=output.init
basenames+=($initfile)
echo ${basenames[*]}

# Rebuild everything in parallel
REBUILD_NEMO=/nobackup/smhid20/users/sm_erimu/NEMO/nemo4/tools/REBUILD/rebuild
time parallel --eta -v "${REBUILD_NEMO} -o {}.nc {}_[0-9]*.nc" ::: ${basenames[*]}

echo "cleanup"
for basename in ${basenames[@]}; do
    mkdir -p tmp
    mv ${basename}_[0-9][0-9]*.nc tmp/.
done
