#!/usr/bin/env python

#     |-----------------------------------------------------------|
#     | This module computes timeseries of maximum, 99 percentile |
#     | and average spurious currents after an HPG errors test.   |
#     | The computation is conducted considering both the entire  |
#     | computational domain and only the localisation area.      |
#     |                                                           |
#     | Author: Diego Bruciaferri                                 |
#     | Date and place: 26-08-2024, Met Office, UK                |
#     |-----------------------------------------------------------|


import os
from os.path import join, isfile, basename, splitext
import glob
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import xarray as xr
from dask.diagnostics import ProgressBar

from lib import compute_masks

# ==============================================================================
# Input files
# ==============================================================================

# Folder path containing HPGE spurious currents velocity files 
MAINdir = '/scratch/dbruciaf/GOSI10p1.0-hpge/'
HPGElst = 'r12_r16-r075-r040-r035_it2-r030'
DOMCFG  = '/data/users/dbruciaf/GOSI10_input_files/MEs_novf/025/test/domain_cfg_MEs_novf_4env_2930_r12_r16-r075-r040-r035_it2-r030.nc'

label = 'MEs'

# Name of the zonal and meridional velocity variables
Uvar = 'uo'
Vvar = 'vo'
# Name of the variable to chunk with dask and size of chunks
chunk_var = 'time_counter'
chunk_size = 1

#cols = ["red","blue","dodgerblue","limegreen"]
cols = ["limegreen"]

# ==============================================================================
# Loading domain geometry
ds_dom = xr.open_dataset(DOMCFG).squeeze()
ds_dom = ds_dom.rename({"nav_lev":"z"})

# Computing land-sea masks
ds_dom = compute_masks(ds_dom, merge=True)

e3t = ds_dom["e3t_0"].squeeze()
e2t = ds_dom["e2t"].squeeze()
e1t = ds_dom["e1t"].squeeze()

e1t = e1t.where(ds_dom.tmask==1)
e2t = e2t.where(ds_dom.tmask==1)
e3t = e3t.where(ds_dom.tmask==1)

e1t = e1t.where(msk==1)
e2t = e2t.where(msk==1)
e3t = e3t.where(msk==1)

cel_vol = e1t * e2t * e3t
dom_vol = cel_vol.sum(skipna=True)

HPGEdir = MAINdir + HPGElst

Ufiles = sorted(glob.glob(HPGEdir+'/*grid-U*.nc'))
Vfiles = sorted(glob.glob(HPGEdir+'/*grid-V*.nc'))

v_max_tot = []
v_99p_tot = []
v_avg_tot  = []

for F in range(len(Ufiles)):

    print(Ufiles[F])

    ds_U = xr.open_dataset(Ufiles[F], chunks={chunk_var:chunk_size})
    ds_V = xr.open_dataset(Vfiles[F], chunks={chunk_var:chunk_size})
    U4   = ds_U[Uvar]
    V4   = ds_V[Vvar]

    # rename some dimensions
    U4 = U4.rename({U4.dims[0]: 't', U4.dims[1]: 'z'})
    V4 = V4.rename({V4.dims[0]: 't', V4.dims[1]: 'z'})

    # interpolating from U,V to T
    U = U4.rolling({'x':2}).mean().fillna(0.)
    V = V4.rolling({'y':2}).mean().fillna(0.) 
    
    vel_t = np.sqrt(np.power(U,2) + np.power(V,2)) 
    vel_t = vel_t.where(ds_dom.tmask==1)

    v_max_tot.extend(vel_t.max(dim=('z','y','x'),skipna=True).values.tolist())
    v_99p_tot.extend(vel_t.quantile(0.99,dim=('z','y','x'),skipna=True).values.tolist())
    v_avg_tot.extend(((cel_vol*vel_t).sum(dim=["x","y","z"], skipna=True) / dom_vol).values.tolist())

# Saving 

ds = xr.Dataset()

ds["max_u_tot"] = xr.DataArray(np.asarray(v_max_tot), dims=('t'))
ds["u_99p_tot"] = xr.DataArray(np.asarray(v_99p_tot), dims=('t'))
ds["avg_u_tot"] = xr.DataArray(np.asarray(v_avg_tot), dims=('t'))

# -------------------------------------------------------------------------------------   
# Writing the max_hpge file

print('WRITING the maximum_hpge.nc FILE')

out_file = "hpge_timeseries.nc"
delayed_obj = ds.to_netcdf(join(HPGEdir,out_file), compute=False)

with ProgressBar():
     results = delayed_obj.compute()

