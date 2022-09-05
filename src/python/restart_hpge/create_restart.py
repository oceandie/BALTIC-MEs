#!/usr/bin/env python

import numpy as np
import xarray as xr
import netCDF4 as nc4
import gsw as gsw
from scipy import interpolate

Tfile  = '/data/users/dbruciaf/mod_dev/dbruciaf-gulf18/vgrid/restart_hpge/sf12-r30_sm5_r8-H20/20140101_restart_T-S_hpge.nc'
Sfile  = '/data/users/dbruciaf/mod_dev/dbruciaf-gulf18/vgrid/restart_hpge/sf12-r30_sm5_r8-H20/20140101_restart_T-S_hpge.nc'
Tname  = 'tn'
Sname  = 'sn'
inpmsh = '/data/users/dbruciaf/mod_dev/files4bdy_gulf18/gulf18.mesh_mask.sf12-r30_sm5_r8-H20.nc'
outmsh = '/data/users/dbruciaf/mod_dev/files4bdy_gulf18/gulf18.mesh_mask.MEs_2env_0.3_0.1_opt.nc'

i1D = 570
j1D = 87

# -----------------------------------------------------------------

# Loading arrays
dsT = xr.open_dataset(Tfile)
pot_tem = dsT[Tname].squeeze().values

dsS = xr.open_dataset(Sfile)
pra_sal = dsS[Sname].squeeze().values

ds_inp = xr.open_dataset(inpmsh)
gdeptI = ds_inp['gdept_0'].squeeze().values
mbkt = ds_inp['mbathy'].squeeze().values

ds_out = xr.open_dataset(outmsh)
gdeptO = ds_out['gdept_0'].squeeze().values

nk = pot_tem.shape[0]
nj = pot_tem.shape[1]
ni = pot_tem.shape[2]

# Extracting the 1D profiles
T1D = pot_tem[:,j1D,i1D]
S1D = pra_sal[:,j1D,i1D]
dep1D = gdeptI[:,j1D,i1D]
wet1D = mbkt[j1D,i1D]

T1D[wet1D:] = T1D[wet1D] # getting rid of nans at depth
S1D[wet1D:] = S1D[wet1D] # getting rid of nans at depth

# Interpolating the 1D profile in the out domain
Tout = np.zeros(shape=gdeptO.shape)
Sout = np.zeros(shape=gdeptO.shape)
fT = interpolate.interp1d(dep1D, T1D, fill_value="extrapolate")
fS = interpolate.interp1d(dep1D, S1D, fill_value="extrapolate")
for j in range(nj):
    for i in range(ni):
        Tout[:,j,i] = fT(gdeptO[:,j,i])
        Sout[:,j,i] = fS(gdeptO[:,j,i])

dsT[Tname].values[0,:,:,:] = Tout[:,:,:]
dsT[Sname].values[0,:,:,:] = Sout[:,:,:]

out_file = './TS_data.nc'
dsT.to_netcdf(out_file)

