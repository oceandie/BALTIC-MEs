#!/usr/bin/env python

import os
import sys
import subprocess
import numpy as np
import xarray as xr
from xnemogcm import open_domain_cfg
from plot_section import mpl_sec_loop
from utils import compute_masks

# ========================================================================
# INPUT PARAMETERS

DOMCFG_zps = '/nobackup/smhid20/users/sm_erimu/NEMO/generate_domcfg/nordic_ref/domain_cfg.nc'
latlonfile = '/nobackup/smhid20/users/sm_erimu/NEMO/generate_domcfg/analysis/lonslats.nc'

# 2. ANALYSIS

sec_lon1 = [13.124825, 13.680375, 14.235925, 14.791475, 15.347025, 15.902575,
            16.458125, 17.013675, 17.569225, 18.124775, 18.569215, 18.84699 ,
            19.29143 , 19.680315, 19.902535, 20.013645]

sec_lon2 = [8.902645, 8.9582  , 9.06931 , 9.124865, 9.18042 ]

sec_lat1 = [54.891598, 55.058262, 55.258258, 55.39159 , 55.458255, 55.291591,
            55.291591, 55.291591, 55.358257, 55.458255, 55.691585, 56.024913,
            56.358241, 56.691569, 57.024897, 57.358225]

sec_lat2 = [58.424874, 58.091546, 57.758218, 57.42489 , 57.091562]

ds = xr.open_dataset(latlonfile)
sec_lat3 = ds.gphit.values.squeeze().tolist()
sec_lon3 = ds.glamt.values.squeeze().tolist()

print(sec_lat2)
print(sec_lon2)
print(sec_lat3)
print(sec_lon3)

# quit()

sec_I_indx_1b_L  = [sec_lon1, sec_lon2, sec_lon3]
sec_J_indx_1b_L  = [sec_lat1, sec_lat2, sec_lat3]

coord_type_1b_L  = "dist"
rbat2_fill_1b_L  = "false"
xlim_1b_L        = "maxmin"
ylim_1b_L        = [0., 260.] #1000.] #3500.] #5900.]
vlevel_1b_L      = 'Z_ps'
xgrid_1b_L       = "false"

# ========================================================================
# Loading domain geometry
# ds_dom  = open_domain_cfg(files=[DOMCFG_zps])
ds_dom = xr.open_dataset(DOMCFG_zps).squeeze()

# for i in ['bathymetry','bathy_meter']:
#     for dim in ['x','y']:
#         ds_dom[i] = ds_dom[i].rename({dim: dim+"_c"})

# Computing masks
ds_dom = compute_masks(ds_dom, merge=True)

tlon2 = ds_dom["glamt"].values
tlat2 = ds_dom["gphit"].values
e3t_3 = ds_dom["e3t_0"].values
e3w_3 = ds_dom["e3w_0"].values
tmsk3 = ds_dom["tmask"].values
bathy = []
#bathy = ds_dom["bathymetry"].values

nk = e3t_3.shape[0]
nj = e3t_3.shape[1]
ni = e3t_3.shape[2]

tlon3 = np.repeat(tlon2[np.newaxis, :, :], nk, axis=0)
tlat3 = np.repeat(tlat2[np.newaxis, :, :], nk, axis=0)

# Computing model levels' depth
tdep3 = np.zeros(shape=(nk,nj,ni))
wdep3 = np.zeros(shape=(nk,nj,ni))
wdep3[0,:,:] = 0.
tdep3[0,:,:] = 0.5 * e3w_3[0,:,:]
for k in range(1, nk):
    wdep3[k,:,:] = wdep3[k-1,:,:] + e3t_3[k-1,:,:]
    tdep3[k,:,:] = tdep3[k-1,:,:] + e3w_3[k,:,:]

proj = []

# PLOTTING VERTICAL DOMAIN

var_strng  = ""
unit_strng = ""
date       = ""
timeres_dm = ""
timestep   = []
PlotType   = ""
var4       = []
hbatt      = []
mbat_ln    = "false"
mbat_fill  = "true"
varlim     = "no"
check      = 'true'
check_val  = 'false'


mpl_sec_loop('AMM15 mesh', '.png', var_strng, unit_strng, date, timeres_dm, timestep, PlotType,
              sec_I_indx_1b_L, sec_J_indx_1b_L, tlon3, tlat3, tdep3, wdep3, tmsk3, var4, proj,
              coord_type_1b_L, vlevel_1b_L, bathy, hbatt, rbat2_fill_1b_L, mbat_ln, mbat_fill,
              xlim_1b_L, ylim_1b_L, varlim, check, check_val, xgrid_1b_L)


