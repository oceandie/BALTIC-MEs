
# Table of Contents

1.  [Compilation](#orgd521709)
    1.  [Compiling `NEMO_4.0.4_hpge` on `Freja` ](#orgccb1c39)
    2.  [NEMO 4.2.2](#org5ea693d)
    3.  [Compiling DOMAINcfg tool](#org19709fb)
2.  [Generate envelopes ](#orge9b746c)
    1.  [Notes](#orge3b0456)
3.  [Create a `domain_cfg.nc`](#orgc806800)
    1.  [visualize `domain_cfg.nc`](#org25fa458)
4.  [HPG error testing](#orgeac42af)
    1.  [(optional) add an initial TS depth-profile (same procedure for 4.2.2)](#org36c4871)
    2.  [run HPGE test](#orgcf45fda)
    3.  [Create `maximum_hpge.nc`](#org7bd234d)
    4.  [HPGE iteration](#org09071ff)


<a id="orgd521709"></a>

# Compilation


<a id="orgccb1c39"></a>

## Compiling `NEMO_4.0.4_hpge` on `Freja` <a id="orgdc35d6e"></a>

On an interactive node

    module purge
    module load buildenv-intel/2023a-eb
    module load netCDF-HDF5/4.9.2-1.12.2-hpc1
    
    cd <REPO_DIR>/src/f90/NEMO_4.0.4_hpge_ovf/;
    ./makenemo -m 'freja' -r HPGTEST -n 'HPGTEST' -j 16


<a id="org5ea693d"></a>

## NEMO 4.2.2

In the fouo nemo4 repo there's a branch available called `hpg_test_422_fixes`
Build it as follows:

    module purge;
    module load buildenv-intel/2023a-eb;
    module load netCDF-HDF5/4.9.2-1.12.2-hpc1;
    ./makenemo -m 'freja' -r NORDIC -n 'NORDIC_HPGTEST_FREJA_422_fixes' -j 16

And use the executable in `nemo4/cfgs/NORDIC_HPGTEST_FREJA_422_fixes/BLD/bin/nemo.exe` for the hpg tests


<a id="org19709fb"></a>

## Compiling DOMAINcfg tool

NEMO 4.2.x domaintool comes with MEs support. From the repo
<https://forge.nemo-ocean.eu/nemo/nemo.git> we use the branch
`423-adding-more-flexibility-to-me-gvcs`

    cd <NEMO_4.2.2_REPO_DIR>/tools
    git checkout 423-adding-more-flexibility-to-me-gvcs
    ./maketools -m 'freja' -n DOMAINcfg


<a id="orge9b746c"></a>

# Generate envelopes <a id="orgb77d3bb"></a>

To generate envelopes:

-   create some dir, for instance here `<REPO_DIR>/src/python/envelopes/<some_dir>`
-   add an existing `bathy_meter.nc` and corresponding `domain_cfg.nc` or `coordinates.nc`
-   create `<inputname>.inp` and set it up
-   link to `<REPO_DIR>/src/python/envelopes/generate_envelopes.py`
-   `python generate_envelopes.py <inputname>.inp`
    probably have to fix missing packages

The result is a new bathymetry file: `bathy_meter.<inputname>.nc`


<a id="orge3b0456"></a>

## Notes

Envelope generation with `generate_envelopes.py` and an input `.inp` file:

-   Stick to the format given in example inp files.
-   Do not change the variable types (float/int).
-   Use only a single value for the threshold, even though the list allows several.


<a id="orgc806800"></a>

# Create a `domain_cfg.nc`

With the new bathymetry file we create a `domain_cfg.nc`, where
additional settings should be given in the `namelist_cfg`.

1.  Create a `coordinates.nc` from a `mesh_mask.nc`. This is only
    for horizontal coordinates so it can be done once using the
    mesh<sub>mask</sub> corresponding to the original ZPS domain:
    
        ncks -v glamt,glamu,glamv,glamf,gphit,gphiu,gphiv,gphif,e1t,e1u,e1v,e1f,e2t,e2u,e2v,e2f,ff_t,ff_f,nav_lon,nav_lat mesh_mask.nc coordinates.nc=

2.  Copy contents of `<REPO_DIR>/src/scripts/create_domaincfg/` to
    some `<DOMAIN_BLD>` dir. To see what I have been running look in
    `<REPO_DIR>/customized_scripts/create_domaincfg`.

3.  Edit `<DOMAIN_BLD>/files/namelist_cfg` to suit your needs. The
    default setup in the repo is for 2 envelopes with 43 and 13
    levels respectively.

4.  On an interactive node,
    
        module purge
        module load buildenv-intel/2023a-eb
        module load netCDF-HDF5/4.9.2-1.12.2-hpc1

5.  Run `./create_domaincfg.sh <outputdir> <bathymetry>`
    this should take less than 2 minutes.


<a id="org25fa458"></a>

## visualize `domain_cfg.nc`

in `<REPO_DIR>/src/python/plot/vcoord/`
edit and run `plot_vlevels_MEs.py` or `plot_vlevels_zps.py`


<a id="orgeac42af"></a>

# HPG error testing


<a id="org36c4871"></a>

## (optional) add an initial TS depth-profile (same procedure for 4.2.2)

1.  add an initial TS depth-profile to
    `<REPO_DIR>/src/f90/NEMO_4.0.4_hpge_ovf/src/OCE/USR/usrdef_istate.F90`
    if necessary.
2.  recompile `NEMO_4.0.4_hpge` ([1.1](#orgdc35d6e))


<a id="orgcf45fda"></a>

## run HPGE test

1.  Copy contents of `<REPO_DIR>/src/scripts/run_hpgtest/` to
    preferred rundir. To see what I have been running look in
    `<REPO_DIR>/customized_scripts/run_hpg_test`.

2.  select initial TS depth-profile in `test_template/namelist_cfg` (`namtsd::nn_tsd_type`)
3.  create rundir with hpge setup and submit
    `./run_hpgetest.sh <testname> <domcfg>`


<a id="org7bd234d"></a>

## Create `maximum_hpge.nc`

-   edit and run `create_2D_hpge_field.py` (in `<REPO_DIR>/src/python/envelopes`)
-   (optional) visualize in the test dir: `ncview maximum_hpge.nc`


<a id="org09071ff"></a>

## HPGE iteration

Not happy with the HPGE? Go back to [2](#orgb77d3bb) and use
 `maximum_hpge.nc` to create a new bathymetry with HPGE aware local
 smoothing (see example `.inp` files). Note that several
 `maximum_hpge.nc` input fields can be used.

Otherwise you're done and you can start running experiments.

