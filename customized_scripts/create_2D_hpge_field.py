#!/usr/bin/env python

#     |------------------------------------------------------------|
#     | This module creates a 2D field of maximum spurious current |
#     | in the vertical and in time after an HPGE test.            |
#     | The resulting file can be used then to optimise the rmax   |
#     | of Multi-Envelope vertical grids.                          |
#     |                                                            |
#     | Author: Diego Bruciaferri                                  |
#     | Date and place: 07-09-2021, Met Office, UK                 |
#     |------------------------------------------------------------|
# Edits made by Erik, SMHI, for the Baltic MEs project/

import os
from os.path import join, isfile, basename, splitext
import glob
import numpy as np
from matplotlib import pyplot as plt
import xarray as xr
from dask.diagnostics import ProgressBar

# ==============================================================================
# Input files
# ==============================================================================

# Folder path containing HPGE spurious currents velocity files
HPGEdir = '/nobackup/smhid20/users/sm_erimu/BALTIC-MEs/hpg_test/emodnet_domcfgs_2envs_ls_1cm_v2_long'
HPGEdir = '/nobackup/smhid20/users/sm_erimu/BALTIC-MEs/hpg_test/emodnet_domcfgs_2envs_v2'

def create_maximum_hpge(add_more=False):

    # List of indexes of the last T-level of each vertical subdomains
    # (Fortran indexening convention)
    num_lev = [43, 56]

    # Name of the zonal and meridional velocity variables
    Uvar = 'uoce_inst'
    Vvar = 'voce_inst'

    # ==============================================================================
    # LOOP

    Ufiles = sorted(glob.glob(HPGEdir+'/*grid_U_*.nc'))
    Vfiles = sorted(glob.glob(HPGEdir+'/*grid_V_*.nc'))

    hpge_lists = dict()

    for F in range(len(Ufiles)):

        print(Ufiles[F])
        print(Vfiles[F])

        ds_U = xr.open_dataset(Ufiles[F])
        U4   = ds_U[Uvar].load()
        ds_V = xr.open_dataset(Vfiles[F])
        V4   = ds_V[Vvar].load()

        # rename some dimensions
        U4 = U4.rename({U4.dims[0]: 't', U4.dims[1]: 'k'})
        V4 = V4.rename({V4.dims[0]: 't', V4.dims[1]: 'k'})

        # interpolating from U,V to T
        U = U4.rolling({'x':2}).mean().fillna(0.)
        V = V4.rolling({'y':2}).mean().fillna(0.)

        hpge = np.sqrt(np.power(U,2) + np.power(V,2))

        if F == 0:
           ni = hpge.data.shape[3]
           nj = hpge.data.shape[2]
           if len(num_lev) > 1:
              max_hpge1 = np.zeros(shape=(nj,ni))
              max_hpge2 = np.zeros(shape=(nj,ni))
              # max_hpge3 = np.zeros(shape=(nj,ni))
              # max_hpge4 = np.zeros(shape=(nj,ni))
              hpge_lists['1'] = []
              hpge_lists['2'] = []
           else:
              max_hpge1 = np.zeros(shape=(nj,ni))

        if len(num_lev) > 1:
            hpge_1 = hpge.isel(k=slice(None, num_lev[0])).max(dim='k')
            hpge_2 = hpge.isel(k=slice(num_lev[0], num_lev[1])).max(dim='k')
            hpge_lists['1'].append(hpge_1)
            hpge_lists['2'].append(hpge_2)

            maxhpge_1 = hpge_1.max(dim='t')
            maxhpge_2 = hpge_2.max(dim='t')
            # maxhpge_3 = hpge.isel(k=slice(num_lev[1], num_lev[2])).max(dim='k').max(dim='t')
            # maxhpge_4 = hpge.isel(k=slice(num_lev[2], num_lev[3])).max(dim='k').max(dim='t')
            max_hpge1 = np.maximum(max_hpge1, maxhpge_1.data)
            max_hpge2 = np.maximum(max_hpge2, maxhpge_2.data)
            # max_hpge3 = np.maximum(max_hpge3, maxhpge_3.data)
            # max_hpge4 = np.maximum(max_hpge4, maxhpge_4.data)
        else:
            hpge_1 = hpge.isel(k=slice(None, num_lev[0])).max(dim='k')
            hpge_lists['1'].append(hpge_1)
            maxhpge_1 = hpge_1.max(dim='t')
            max_hpge1 = np.maximum(max_hpge1, maxhpge_1.data)


    # Saving
    ds_hpge = xr.Dataset()
    if len(num_lev) > 1:
       ds_hpge["max_hpge_1"] = xr.DataArray(max_hpge1, dims=('y','x'))
       ds_hpge["max_hpge_2"] = xr.DataArray(max_hpge2, dims=('y','x'))
       # ds_hpge["max_hpge_3"] = xr.DataArray(max_hpge3, dims=('y','x'))
       # ds_hpge["max_hpge_4"] = xr.DataArray(max_hpge4, dims=('y','x'))
    else:
       ds_hpge["max_hpge_1"] = xr.DataArray(max_hpge1, dims=('y','x'))

    if add_more:
        ds_hpge["hpge_1"] = xr.concat(hpge_lists['1'], dim='t')
        if len(num_lev) > 1:
            ds_hpge["hpge_2"] = xr.concat(hpge_lists['2'], dim='t')

    # -------------------------------------------------------------------------------------
    # Writing the max_hpge file

    if add_more:
        out_file = "maximum_hpge_extended.nc"
    else:
        out_file = "maximum_hpge.nc"

    print(f'Loading into memory')
    ds_hpge.load()
    print(f'WRITING to {out_file}')
    delayed_obj = ds_hpge.to_netcdf(join(HPGEdir, out_file), compute=False)

    with ProgressBar():
         results = delayed_obj.compute()

if __name__=='__main__':
    for flag in [False, True]:
        create_maximum_hpge(flag)
