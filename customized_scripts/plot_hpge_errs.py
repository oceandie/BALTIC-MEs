import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
plt.close('all')

reference_velocities='/home/sm_erimu/Projects/analysis/data/vel_max_col_mn.pkl'
with open(reference_velocities, 'rb') as fle:
    data = pd.read_pickle(fle)

vel_ref = data['vel_max_col_mn'].sel(src='hindcast ZPS')\
                                .drop_vars(['x','y'])\
                                .fillna(np.inf)
vel_ref = xr.where(vel_ref == 0.0, np.inf, vel_ref)


# v2 HPG error test with NEMO 4.2 vec / bilapl / triad
data_dict = {
    # 'hpge v0' : xr.open_dataset('/nobackup/smhid20/users/sm_erimu/BALTIC-MEs/hpg_test/'
    #                             'emodnet_domcfgs_2envs_v0/'
    #                             'maximum_hpge_extended.nc'),
    # 'hpge 1cm v1' : xr.open_dataset('/nobackup/smhid20/users/sm_erimu/BALTIC-MEs/hpg_test/'
    #                                 'emodnet_domcfgs_2envs_ls_1cm/'
    #                                 'maximum_hpge_extended.nc'),
    # 'hpge 1cm v2' : xr.open_dataset('/nobackup/smhid20/users/sm_erimu/BALTIC-MEs/hpg_test/'
    #                                 'emodnet_domcfgs_2envs_ls_1cm_v2/'
    #                                 'maximum_hpge_extended.nc'),
    # 'hpge 1cm v3' : xr.open_dataset('/nobackup/smhid20/users/sm_erimu/BALTIC-MEs/hpg_test/'
    #                                 'emodnet_domcfgs_2envs_ls_1cm_v3/'
    #                                 'maximum_hpge_extended.nc'),
    # 'hpge v0 r007-007' : xr.open_dataset('/nobackup/smhid20/users/sm_erimu/BALTIC-MEs/hpg_test/'
    #                                      'emodnet_domcfgs_2envs_r007-007/'
    #                                      'maximum_hpge_extended.nc'),
    # 'hpge v0 r007-007 v1' : xr.open_dataset('/nobackup/smhid20/users/sm_erimu/BALTIC-MEs/hpg_test/'
    #                                         'emodnet_domcfgs_2envs_r007-007_ls_1cm/'
    #                                         'maximum_hpge_extended.nc'),
    # 'hpge v0 r007-007 v2' : xr.open_dataset('/nobackup/smhid20/users/sm_erimu/BALTIC-MEs/hpg_test/'
    #                                         'emodnet_domcfgs_2envs_r007-007_ls_1cm_v2/'
    #                                         'maximum_hpge_extended.nc'),
    
    'hpge r007-007 newtool' : xr.open_dataset('/nobackup/smhid20/users/sm_erimu/BALTIC-MEs/hpg_test/'
                                              'emodnet_domcfgs_2envs_r007-007_newtool/'
                                              'maximum_hpge_extended.nc'),
    'hpge r007-007 v1 newtool' : xr.open_dataset('/nobackup/smhid20/users/sm_erimu/BALTIC-MEs/hpg_test/'
                                                 'emodnet_domcfgs_2envs_r007-007_ls_1cm_newtool/'
                                                 'maximum_hpge_extended.nc'),
    'hpge r007-007 v2 newtool' : xr.open_dataset('/nobackup/smhid20/users/sm_erimu/BALTIC-MEs/hpg_test/'
                                                 'emodnet_domcfgs_2envs_r007-007_ls_1cm_v2_newtool/'
                                                 'maximum_hpge_extended.nc'),
    'hpge r007-007 v3 newtool' : xr.open_dataset('/nobackup/smhid20/users/sm_erimu/BALTIC-MEs/hpg_test/'
                                                 'emodnet_domcfgs_2envs_r007-007_ls_1cm_v3_newtool/'
                                                 'maximum_hpge_extended.nc'),

    
    # 'hpge v0 r007-007 v2 newtool' : xr.open_dataset('/nobackup/smhid20/users/sm_erimu/BALTIC-MEs/hpg_test/'
    #                                         'emodnet_domcfgs_2envs_r007-007_ls_1cm_v2_newtool/'
    #                                         'maximum_hpge_extended.nc'),
    # 'hpge v0 r007-007 v3' : xr.open_dataset('/nobackup/smhid20/users/sm_erimu/BALTIC-MEs/hpg_test/'
    #                                         'emodnet_domcfgs_2envs_r007-007_ls_1cm_v3/'
    #                                         'maximum_hpge_extended.nc'),
    # 'hpge v0 r005-007' : xr.open_dataset('/nobackup/smhid20/users/sm_erimu/BALTIC-MEs/hpg_test/'
    #                                      'emodnet_domcfgs_2envs_r005-007/'
    #                                      'maximum_hpge_extended.nc'),
    # 'hpge v0 r005-007 v2' : xr.open_dataset('/nobackup/smhid20/users/sm_erimu/BALTIC-MEs/hpg_test/'
    #                                         'emodnet_domcfgs_2envs_r005-007_ls_1cm_v2/'
    #                                         'maximum_hpge_extended.nc'),

    # 'hpge 1cm v4' : xr.open_dataset('/nobackup/smhid20/users/sm_erimu/BALTIC-MEs/hpg_test/'
    #                                 'emodnet_domcfgs_2envs_ls_1cm_v4/'
    #                                 'maximum_hpge_extended.nc'),
    # 'hpge 1cm v4 alt1' : xr.open_dataset('/nobackup/smhid20/users/sm_erimu/BALTIC-MEs/hpg_test/'
    #                                      'emodnet_domcfgs_2envs_ls_1cm_v4_alt1/'
    #                                      'maximum_hpge_extended.nc'),
    # 'hpge 1cm v4 alt2' : xr.open_dataset('/nobackup/smhid20/users/sm_erimu/BALTIC-MEs/hpg_test/'
    #                                      'emodnet_domcfgs_2envs_ls_1cm_v4_alt2/'
    #                                     'maximum_hpge_extended.nc'),
    # 'hpge 1cm v5' : xr.open_dataset('/nobackup/smhid20/users/sm_erimu/BALTIC-MEs/hpg_test/'
    #                                 'emodnet_domcfgs_2envs_ls_1cm_v5/'
    #                                 'maximum_hpge_extended.nc'),
    # 'hpge 1cm v6' : xr.open_dataset('/nobackup/smhid20/users/sm_erimu/BALTIC-MEs/hpg_test/'
    #                                 'emodnet_domcfgs_2envs_ls_1cm_v5/'
    #                                 'maximum_hpge_extended.nc'),
    # 'hpge 1.5cm v1' : xr.open_dataset('/nobackup/smhid20/users/sm_erimu/BALTIC-MEs/hpg_test/'
    #                                   'emodnet_domcfgs_2envs_ls_1p5cm/'
    #                                   'maximum_hpge_extended.nc'),
    # 'hpge 1.5cm v2' : xr.open_dataset('/nobackup/smhid20/users/sm_erimu/BALTIC-MEs/hpg_test/'
    #                                   'emodnet_domcfgs_2envs_ls_1p5cm_v2/'
#                                   'maximum_hpge_extended.nc'),
}

sitesfile = '/home/sm_erimu/Projects/NEMO/INPUT_FILES_NEMO4.2.0/grid_cfsites.nc'
dsites = xr.open_dataset(sitesfile, decode_times=False)
by15=dsites.site[9]
by15lat = float(by15.lat)
by15lon = float(by15.lon)

domainfile = '/home/sm_erimu/Projects/BALTIC-MEs/hpg_test/emodnet_domcfgs_2envs_v0/domain_cfg.nc'
ddomain = xr.open_dataset(domainfile)
lats = ddomain.nav_lat[:,0].values
lons = ddomain.nav_lon[0,:].values

bathyfile='/home/sm_erimu/Projects/BALTIC-MEs/src/python/envelopes/nordic_emodnet_domcfg/bathy_meter_EMODnet.v1.0.nordic_mes_2env_r01-007_localsmoothing_v1.nc'
dbath = xr.open_dataset(bathyfile)





def coords_to_inds(lat, lon, lats, lons):
    """
        convert lat,lon-coordinates to grid indices xind,yind
    """
    lonmin = np.amin(lons)
    lonmax = np.amax(lons)
    latmin = np.amin(lats)
    latmax = np.amax(lats)

    lon_in_range = bool((lon >= lonmin) and (lon <= lonmax))
    lat_in_range = bool((lat >= latmin) and (lat <= latmax))

    if (lon_in_range) and (lat_in_range):
        xid = np.argmin(np.abs(lons - lon))
        yid = np.argmin(np.abs(lats - lat))
    else:
        xid = np.nan
        yid = np.nan

    return xid, yid

xid, yid = coords_to_inds(by15lat, by15lon, lats, lons)

plt.close('all')

maxfield='max_hpge_1'#, 'max_hpge_2']
hpgfield='hpge_1'#, 'hpge_2']

numrows = int(np.ceil(len(data_dict)/2) + 1)

for mode in ['normal', 'relative']:


    plt.figure()#figsize=(15,20))
    #plt.subplot(numrows,2,1)

    for key, item in data_dict.items():
        if mode == 'normal':
            field = item[hpgfield].where(np.isnan(vel_ref) == False)
        elif mode == 'relative':
            field = item[hpgfield].where(np.isnan(vel_ref) == False) / vel_ref
            
        plt.plot(field.max(dim=['x','y']), label=key)

    ax = plt.gca()
    ax.set_xlabel('days')
    ax.set_title(f'{hpgfield}')
    ax.legend()

    # plt.subplot(numrows,2,2)

    # for key, item in data_dict.items():
    #     plt.hist(item[hpgfield].data.flatten(), bins=1000, label=key, alpha=0.8)

    # ax = plt.gca()
    # ax.set_title('hpge_1 distribution tail')
    # vmax=0.03 if mode == 'normal' else 0.75
    # ax.set_xlim([0.015,vmax])
    # ax.set_ylim([0.,5e2])

    # vmin=0.0

    # for i, (key, item) in enumerate(data_dict.items()):
    #     plt.subplot(numrows,2,3+i)
        
    #     if mode == 'normal':
    #         field = item[maxfield].where(np.isnan(vel_ref) == False)
    #         label = 'error velocity (m/s)'
    #     elif mode == 'relative':
    #         field = item[maxfield].where(np.isnan(vel_ref) == False) / vel_ref
    #         label = 'relative error velocity'
            
    #     field.plot.imshow(vmin=vmin, vmax=vmax, cbar_kwargs={'label' : label})            
    #     plt.plot(xid, yid,'r.-', markersize=15)
    #     plt.annotate(' by15', [xid,yid], color='r')
    #     ax = plt.gca()
    #     ax.set_title(key)
    #     # ax.set_xlim([300,500])
    #     # ax.set_ylim([150,350])

    plt.tight_layout()
    figname = f'hpge_errs_{hpgfield}.png' if mode == 'normal' \
        else f'hpge_errs_{hpgfield}_{mode}.png'
    plt.savefig(figname)

    plt.figure(figsize=(20,32))

    dirs = [
        'emodnet_domcfgs_2envs_r007-007',
        'emodnet_domcfgs_2envs_r007-007_ls_1cm',
        'emodnet_domcfgs_2envs_r007-007_ls_1cm_v2',
        'emodnet_domcfgs_2envs_r007-007_ls_1cm_v3',
        'emodnet_domcfgs_2envs_r007-007_newtool',
        'emodnet_domcfgs_2envs_r007-007_ls_1cm_newtool',
        'emodnet_domcfgs_2envs_r007-007_ls_1cm_v2_newtool',
        'emodnet_domcfgs_2envs_r007-007_ls_1cm_v3_newtool',
        # 'emodnet_domcfgs_2envs_ls_1cm/',
        # 'emodnet_domcfgs_2envs_ls_1cm_v2/',
        # 'emodnet_domcfgs_2envs_ls_1cm_v3/',
        # 'emodnet_domcfgs_2envs_ls_1cm_v4/',
        #     'emodnet_domcfgs_2envs_ls_1cm_v4_alt1/',
        #     'emodnet_domcfgs_2envs_ls_1cm_v4_alt2/',
        #     'emodnet_domcfgs_2envs_ls_1cm_v5/',
        #     'emodnet_domcfgs_2envs_ls_1cm_v6/',
        # 'emodnet_domcfgs_2envs_ls_1p5cm/',
        # 'emodnet_domcfgs_2envs_ls_1p5cm_v2/',
    ]

    k = 0
    for i, dirnm in enumerate(dirs):
        k += 1
        plt.subplot(len(dirs),2,k)
        try:
            ds = xr.open_dataset('/home/sm_erimu/Projects/BALTIC-MEs/hpg_test'
                                 f'/{dirnm}/maximum_hpge.nc')
        except Exception as e:
            print(e)
            ds = xr.open_dataset('/home/sm_erimu/Projects/BALTIC-MEs/hpg_test'
                                 f'/{dirnm}/maximum_hpge_extended.nc')

        if mode == 'normal':
            field = ds[maxfield].where(np.isnan(vel_ref) == False)
            label = 'error velocity (m/s)'
            add = ''
        elif mode == 'relative':
            field = ds[maxfield].where(np.isnan(vel_ref) == False) / vel_ref
            label = 'relative error velocity'
            add = '(relative)'
            
        field.plot.imshow(vmin=vmin, vmax=vmax,cbar_kwargs={'label' : label})
        plt.plot(xid, yid,'r.-', markersize=15)
        plt.annotate(' by15', [xid,yid], color='r')
        maxval = float(field.max().data)
        print(f'{dirnm}, {maxfield}: {maxval}')
        ax = plt.gca()
        ax.set_title(f'{dirnm}, max={maxval:.8f} {add}', fontsize=10)

        k += 1
        plt.subplot(len(dirs),2,k)
        field.plot.imshow(vmin=vmin, vmax=vmax,cbar_kwargs={'label' : label})
        plt.plot(xid, yid,'r.-', markersize=15)
        plt.annotate(' by15', [xid,yid], color='r')
        ddomain.bathymetry[0,:,:].plot.contour(levels=20, colors='k', alpha=0.5)
        ax = plt.gca()
        ax.set_title(f'{dirnm}', fontsize=10)
        ax.set_xlim([300,500])
        ax.set_ylim([150,350])

    plt.tight_layout()
    figname = f'hpg_iteration_{maxfield}.png' if mode == 'normal' \
        else f'hpg_iteration_{maxfield}_{mode}.png'
    plt.savefig(figname)
