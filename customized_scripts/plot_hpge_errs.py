import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
plt.close('all')

# v2 HPG error test with NEMO 4.2 vec / bilapl / triad
data_dict = {'hpge_v0' : xr.open_dataset('/nobackup/smhid20/users/sm_erimu/BALTIC-MEs/hpg_test/'
                                         'emodnet_domcfgs_2envs_v0/'
                                         'maximum_hpge_extended.nc'),
             'hpge 1cm v1' : xr.open_dataset('/nobackup/smhid20/users/sm_erimu/BALTIC-MEs/hpg_test/'
                                             'emodnet_domcfgs_2envs_ls_1cm/'
                                             'maximum_hpge_extended.nc'),
             'hpge 1cm v2 long' : xr.open_dataset('/nobackup/smhid20/users/sm_erimu/BALTIC-MEs/hpg_test/'
                                                  'emodnet_domcfgs_2envs_ls_1cm_v2_long/'
                                                  'maximum_hpge_extended.nc'),
             'old v2' : xr.open_dataset('/nobackup/smhid20/users/sm_erimu/BALTIC-MEs/hpg_test/'
                                        'emodnet_domcfgs_2envs_v2/'
                                        'maximum_hpge_extended.nc'),
             # 'hpge_v2_2' : xr.open_dataset('/nobackup/smhid20/users/sm_erimu/BALTIC-MEs/hpg_test/'
             #                               'emodnet_domcfgs_2envs_v2_2/'
             #                               'maximum_hpge_extended.nc'),
             }

sitesfile = '/home/sm_erimu/Projects/NEMO/INPUT_FILES_NEMO4.2.0/grid_cfsites.nc'
dsites = xr.open_dataset(sitesfile, decode_times=False)
by15=dsites.site[9]
by15lat = float(by15.lat)
by15lon = float(by15.lon)

domainfile = '/home/sm_erimu/Projects/BALTIC-MEs/hpg_test/emodnet_domcfgs_2envs_v1/domain_cfg.nc'
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
maxfields=['max_hpge_1', 'max_hpge_2']
hpgfields=['hpge_1', 'hpge_2']

for maxfield, hpgfield in zip(maxfields, hpgfields):

    plt.figure(figsize=(15,13))
    plt.subplot(3,2,1)
    
    for key, item in data_dict.items():
        plt.plot(item[hpgfield].max(dim=['x','y']), label=key)

    ax = plt.gca()
    ax.set_xlabel('days')
    ax.set_title(f'{hpgfield}')
    ax.legend()

    plt.subplot(3,2,2)

    for key, item in data_dict.items():
        plt.hist(item[hpgfield].data.flatten(), bins=1000, label=key, alpha=0.5)

    ax = plt.gca()
    ax.set_title('hpge_1 distribution tail')
    vmax=0.03
    ax.set_xlim([0.015,vmax])
    ax.set_ylim([0.,5e2])

    vmin=0.01

    for i, (key, item) in enumerate(data_dict.items()):
        plt.subplot(3,2,3+i)
        item[maxfield].plot.imshow(vmin=vmin, vmax=vmax)
        plt.plot(xid, yid,'r.-', markersize=15)
        plt.annotate(' by15', [xid,yid], color='r')
        ax = plt.gca()
        ax.set_title(key)
        # ax.set_xlim([300,500])
        # ax.set_ylim([150,350])

    plt.tight_layout()
    plt.savefig(f'hpge_errs_{hpgfield}.png')

    plt.figure(figsize=(20,17))

    dirs = ['emodnet_domcfgs_2envs_v0',
            'emodnet_domcfgs_2envs_ls_1cm/',
            'emodnet_domcfgs_2envs_ls_1cm_v2_long/',
            'emodnet_domcfgs_2envs_v2/'
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


        ds[maxfield].plot.imshow(vmin=vmin, vmax=vmax)
        plt.plot(xid, yid,'r.-', markersize=15)
        plt.annotate(' by15', [xid,yid], color='r')
        maxval = float(ds[maxfield].max().data)
        print(f'{dirnm}, {maxfield}: {maxval}')
        ax = plt.gca()
        ax.set_title(f'{dirnm}, max={maxval:.8f}', fontsize=10)

        k += 1
        plt.subplot(len(dirs),2,k)
        ds[maxfield].plot.imshow(vmin=vmin, vmax=vmax)
        plt.plot(xid, yid,'r.-', markersize=15)
        plt.annotate(' by15', [xid,yid], color='r')
        ddomain.bathymetry[0,:,:].plot.contour(levels=20, colors='k', alpha=0.5)
        ax = plt.gca()
        ax.set_title(f'{dirnm}', fontsize=10)
        ax.set_xlim([300,500])
        ax.set_ylim([150,350])

    plt.tight_layout()
    plt.savefig(f'hpg_iteration_{maxfield}.png')    
