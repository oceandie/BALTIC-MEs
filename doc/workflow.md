
# Table of Contents

1.  [Compilation](#org5518935)
    1.  [Compiling `NEMO_4.0.4_hpge` on `Freja` ](#orgaf42b67)
    2.  [Compiling MEs adapted DOMAINcfg tool](#orgbacd898)
    3.  [NEMO 4.2.2](#org7d1a38d)
2.  [Generate envelopes ](#org99f65e4)
    1.  [Notes](#org97b901c)
3.  [Create a `domain_cfg.nc`](#org1424c6a)
    1.  [`namelist_ref::namzgr_mes` for the DOMAINcfg tool ](#org6a98d0c)
    2.  [visualize `domain_cfg.nc`](#orgb3addd0)
4.  [HPG error testing](#orgc7908d2)
    1.  [(optional) add an initial TS depth-profile (same procedure for 4.2.2)](#org8a39ee5)
    2.  [run HPGE test](#org27bce9b)
    3.  [Create `maximum_hpge.nc`](#org884f0ca)
    4.  [HPGE iteration](#org245f5b4)


<a id="org5518935"></a>

# Compilation


<a id="orgaf42b67"></a>

## Compiling `NEMO_4.0.4_hpge` on `Freja` <a id="orgc578959"></a>

On an interactive node

    module purge
    module load buildenv-intel/2023a-eb
    module load netCDF-HDF5/4.9.2-1.12.2-hpc1
    
    cd <REPO_DIR>/src/f90/NEMO_4.0.4_hpge_ovf/;
    ./makenemo -m 'freja' -r HPGTEST -n 'HPGTEST' -j 16


<a id="orgbacd898"></a>

## Compiling MEs adapted DOMAINcfg tool

    cd <REPO_DIR>/src/f90/NEMO_4.0.4_hpge_ovf/tools
    ./maketools -m 'freja' -n DOMAINcfg


<a id="org7d1a38d"></a>

## NEMO 4.2.2

In the fouo nemo4 repo there's a branch available called `hpg_test_422_fixes`
Build it as follows:

    module purge;
    module load buildenv-intel/2023a-eb;
    module load netCDF-HDF5/4.9.2-1.12.2-hpc1;
    ./makenemo -m 'freja' -r NORDIC -n 'NORDIC_HPGTEST_FREJA_422_fixes' -j 16

And use the executable in `nemo4/cfgs/NORDIC_HPGTEST_FREJA_422_fixes/BLD/bin/nemo.exe` for the hpg tests


<a id="org99f65e4"></a>

# Generate envelopes <a id="org1ba855d"></a>

To generate envelopes:

-   create some dir, for instance here `<REPO_DIR>/src/python/envelopes/<some_dir>`
-   add an existing `bathy_meter.nc` and corresponding `domain_cfg.nc` or `coordinates.nc`
-   create `<inputname>.inp` and set it up
-   link to `<REPO_DIR>/src/python/envelopes/generate_envelopes.py`
-   `python generate_envelopes.py <inputname>.inp`
    probably have to fix missing packages

The result is a new bathymetry file: `bathy_meter.<inputname>.nc`


<a id="org97b901c"></a>

## Notes

Envelope generation with `generate_envelopes.py` and an input `.inp` file:

-   Stick to the format given in example inp files.
-   Do not change the variable types (float/int).
-   Use only a single value for the threshold, even though the list allows several.


<a id="org1424c6a"></a>

# Create a `domain_cfg.nc`

With the new bathymetry file we create a `domain_cfg.nc`, where
additional settings should be given in the `namelist_cfg`.

1.  Copy contents of `<REPO_DIR>/src/scripts/create_domaincfg/` to some `<DOMAIN_BLD>` dir.
2.  Edit `<DOMAIN_BLD>/files/namelist_cfg` to suit your needs (see [3.1](#orgf09c573)).
3.  On an interactive node,
    
        module purge
        module load buildenv-intel/2023a-eb
        module load netCDF-HDF5/4.9.2-1.12.2-hpc1
4.  Run `./create_domaincfg.sh <outputdir> <bathymetry>`
    this should take less than 2 minutes.


<a id="org6a98d0c"></a>

## `namelist_ref::namzgr_mes` for the DOMAINcfg tool <a id="orgf09c573"></a>

In this example for NEMO-NORDIC we use two envelopes with 43 and 13 layers respectively.

    !-----------------------------------------------------------------------
    &namzgr_mes    !   MEs-coordinate
    !-----------------------------------------------------------------------
       ln_envl     =   .TRUE. , .TRUE. , .FALSE. , .FALSE., .FALSE.  ! (T/F) If the envelope is used
       nn_strt     =     2    ,   1    ,    1   ,   1    ,   1     ! Stretch. funct.: Madec 1996 (0) or
                                                                   ! Song & Haidvogel 1994 (1) or
                                                                   ! Siddorn & Furner 2012 (2)
       nn_slev     =     43    ,   13   ,   20   ,   15   ,    0   ! number of s-lev between env(n-1)
                                                                   ! and env(n)
       rn_e_hc     =     25.0  ,   0.0 ,   0.0  ,   0.0  ,   0.0   ! critical depth for transition to
                                                                   ! stretch. coord.
       rn_e_th     =     1.5  ,    2.5 ,   2.4  ,   0.0  ,   0.0   ! surf. control param.:
                                                                   ! SH94 or MD96: 0<=th<=20
                                                                   ! SF12: thickness surf. cell
       rn_e_bb     =     -0.3 ,    0.5 ,   0.85 ,   0.0  ,   0.0   ! bot. control param.:
                                                                   ! SH94 or MD96: 0<=bb<=1
                                                                   ! SF12: offset for calculating Zb
       rn_e_al     =     2.0  ,     0.0 ,   0.0  ,   0.0  ,   0.0   ! alpha stretching param with SF12
       rn_e_ba     =     0.024  ,    0.0 ,   0.0  ,   0.0  ,   0.0   ! SF12 bathymetry scaling factor for
                                                                   ! calculating Zb
       rn_bot_min  = 9.0       ! minimum depth of the ocean bottom (>0) (m)
       rn_bot_max  = 700.0     ! maximum depth of the ocean bottom (= ocean depth) (>0) (m)
    
       ln_pst_mes   = .false.
       ln_pst_l2g   = .false.
       rn_e3pst_min = 20.
       rn_e3pst_rat = 0.1
    /


<a id="orgb3addd0"></a>

## visualize `domain_cfg.nc`

in `<REPO_DIR>/src/python/plot/vcoord/`
edit and run `plot_vlevels_MEs.py` or `plot_vlevels_zps.py`


<a id="orgc7908d2"></a>

# HPG error testing


<a id="org8a39ee5"></a>

## (optional) add an initial TS depth-profile (same procedure for 4.2.2)

1.  add an initial TS depth-profile to
    `<REPO_DIR>/src/f90/NEMO_4.0.4_hpge_ovf/src/OCE/USR/usrdef_istate.F90`
    if necessary.
2.  recompile `NEMO_4.0.4_hpge` ([1.1](#orgc578959))


<a id="org27bce9b"></a>

## run HPGE test

1.  copy contents of `<REPO_DIR>/src/scripts/run_hpgtest/` to preferred rundir
2.  select initial TS depth-profile in `test_template/namelist_cfg` (`namtsd::nn_tsd_type`)
3.  create rundir with hpge setup and submit
    `./run_hpgetest.sh <testname> <domcfg>`


<a id="org884f0ca"></a>

## Create `maximum_hpge.nc`

-   edit and run `create_2D_hpge_field.py` (in `<REPO_DIR>/src/python/envelopes`)
-   (optional) visualize in the test dir: `ncview maximum_hpge.nc`


<a id="org245f5b4"></a>

## HPGE iteration

Not happy with the HPGE? Go back to [2](#org1ba855d) and use
 `maximum_hpge.nc` to create a new bathymetry with HPGE aware local
 smoothing (see example `.inp` files). Note that several
 `maximum_hpge.nc` input fields can be used.

Otherwise you're done and you can start running experiments.

