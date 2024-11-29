import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from dask.diagnostics import ProgressBar
from importlib import reload
import lib
import os
reload(lib)
from lib import compute_masks

reference_velocities='/home/sm_erimu/Projects/analysis/data/vel_max_col_mn.pkl'
with open(reference_velocities, 'rb') as fle:
    data = pd.read_pickle(fle)

vel_ref = data['vel_max_col_mn'].sel(src='hindcast ZPS')\
                                .drop_vars(['x','y'])\
                                .fillna(np.inf)
vel_ref = xr.where(vel_ref == 0.0, np.inf, vel_ref)

plt.close('all')
dir_array = [
    'emodnet_domcfgs_2envs_v0',
    'emodnet_domcfgs_2envs_ls_1cm',
    'emodnet_domcfgs_2envs_ls_1cm_v2',
    # 'emodnet_domcfgs_2envs_ls_1cm_v3',
    'emodnet_domcfgs_2envs_ls_1cm_v3_alt1',
    'emodnet_domcfgs_2envs_r007-007',
    'emodnet_domcfgs_2envs_r007-007_ls_1cm',
    'emodnet_domcfgs_2envs_r007-007_ls_1cm_v2',
    'emodnet_domcfgs_2envs_r007-007_ls_1cm_v3',
    # 'emodnet_domcfgs_2envs_r007-007_ls_05cm',
    'emodnet_domcfgs_2envs_r007-007_newtool',
    'emodnet_domcfgs_2envs_r007-007_ls_1cm_newtool',
    'emodnet_domcfgs_2envs_r007-007_ls_1cm_v2_newtool',
    'emodnet_domcfgs_2envs_r007-007_ls_1cm_v3_newtool',
    #'emodnet_domcfgs_2envs_r007-007_ls_1cm_v2_newtool',
    # 'emodnet_domcfgs_2envs_r005-007',
    # 'emodnet_domcfgs_2envs_r005-007_ls_1cm',
    # 'emodnet_domcfgs_2envs_r005-007_ls_1cm_v2',
    # 'emodnet_domcfgs_2envs_r003-007',
    # 'emodnet_domcfgs_2envs_r004-004',
    # 'emodnet_domcfgs_2envs_ls_1cm_v4',
    # 'emodnet_domcfgs_2envs_ls_1cm_v4_alt1',
    # 'emodnet_domcfgs_2envs_ls_1cm_v4_alt2',
    # 'emodnet_domcfgs_2envs_ls_1cm_v5',
    # 'emodnet_domcfgs_2envs_ls_1cm_v6',
    # 'emodnet_domcfgs_2envs_ls_1p5cm',
    # 'emodnet_domcfgs_2envs_ls_1p5cm_v2'
]

mode == 'normal'

for data_dir in dir_array:
    if (
            not data_dir == 'emodnet_domcfgs_2envs_r007-007_ls_1cm_v2_newtool'
    ):
        continue

    print(f'{data_dir}')
    domcfg=f'{data_dir}/domain_cfg.nc'
    
    newtool_mode=True if 'newtool' in data_dir else False
    
    if newtool_mode:
        ds_dom = xr.open_dataset(domcfg).squeeze().rename_dims({'nav_lev':'z'})
    else:
        ds_dom = xr.open_dataset(domcfg).squeeze()
        
    ds_dom = compute_masks(ds_dom, merge=True)

    if newtool_mode:
        tmask = ds_dom.tmask\
                      .drop_vars(['nav_lev'])\
                      .assign_coords({'z':ds_dom.nav_lev.data})
    else:
        tmask = ds_dom.tmask

    if newtool_mode:
        e3t = ds_dom["e3t_0"].squeeze()\
                             .drop_vars(['nav_lev'])\
                             .assign_coords({'z':ds_dom.nav_lev.data})
    else:
        e3t = ds_dom["e3t_0"].squeeze()

    e2t = ds_dom["e2t"].squeeze()
    e1t = ds_dom["e1t"].squeeze()

    e1t = e1t.where(tmask==1)
    e2t = e2t.where(tmask==1)
    e3t = e3t.where(tmask==1)

    e1t = e1t.where(tmask==1)
    e2t = e2t.where(tmask==1)
    e3t = e3t.where(tmask==1)

    cel_vol = e1t * e2t * e3t
    dom_vol = cel_vol.sum(skipna=True)

    chunk_size = 1
    da_u = xr.open_mfdataset(f'{data_dir}/*grid_U*.nc', parallel=True,
                             chunks={'time_counter':chunk_size})\
             .uoce_inst.rename({'depthu' : 'z'})
    da_v = xr.open_mfdataset(f'{data_dir}/*grid_V*.nc', parallel=True,
                             chunks={'time_counter':chunk_size})\
             .voce_inst.rename({'depthv' : 'z'})

    vel = np.sqrt(np.square(da_u) + np.square(da_v))    

    vel_t = vel.max(axis=(1,2,3))

    print('load velocities')
    with ProgressBar():
        vel_t.load()

    u_roll = da_u.rolling({'x':2}).mean().fillna(0.)
    v_roll = da_v.rolling({'y':2}).mean().fillna(0.)

    vel_roll = np.sqrt(np.square(u_roll) + np.square(v_roll))

    # with ProgressBar():
    #     vel_roll_t = vel_roll.max(axis=(1,2,3)).load()

    # # mask-out the North Sea, Kattegat and Sound
    # tmask[...,:310]=0
    # tmask[...,:10,:]=0
    vel_roll_masked = vel_roll.where(tmask==1)

    vel_roll_masked = vel_roll_masked / vel_ref
    
    print('compute rolling velocities')
    with ProgressBar():
        vel_roll_masked_t = vel_roll_masked.max(axis=(1,2,3)).load()

    print('compute 99 percent quantile')
    with ProgressBar():
        vel_roll_masked_q = \
            vel_roll_masked.quantile(0.99,dim=('z','y','x'),skipna=True).load()

    print('compute mean')
    with ProgressBar():
        vel_roll_masked_mn = \
            ((cel_vol*vel_roll_masked).sum(dim=('z','y','x'),skipna=True) / dom_vol).load()

    plt.figure()
    # plt.plot(vel_t.data, label='max velocity')
    add = 'relative ' if mode == 'relative' else ''
    plt.plot(vel_roll_masked_t.data, label=f'max {add}velocity')
    plt.plot(vel_roll_masked_q.data, label=f'99 percent quantile')
    plt.plot(vel_roll_masked_mn.data, label=f'mean {add}velocity')
    plt.legend()
    ylab = 'relative error' if mode == 'relative' else 'velocity (m/s)'
    plt.gca().set_ylabel(ylab)
    # plt.gca().set_ylim([0,0.035])

    title = f'hpg_err_{data_dir}_relative' if mode == 'relative' else f'hpg_err_{data_dir}'
    fname = f'{title}.png'
    plt.suptitle(fname)
    plt.grid()
    print(fname)
    plt.savefig(fname)
    
    # plt.pause(1)
    # m_fname=f'{fname[:-4]}_masked.png'
    # plt.savefig(m_fname)
    # plt.close('all')
    # m_fname2d=f'{fname[:-4]}_masked_2d.png'
    # vel_roll_masked[0,].max(axis=0).plot()
    # plt.savefig(m_fname2d)
    # breakpoint()
    # raise Exception('exceptional')

pngfiles = [f'hpg_err_{dir}.png' for dir in dir_array]
cmd = f"montage {' '.join(pngfiles)} -tile 4x6 -geometry +0+0 all_timeseries.png"
print(cmd)
os.system(cmd)
