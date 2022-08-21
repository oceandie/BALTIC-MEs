# BALTIC-MEs
Collaboration project between the Swedish Meteorological and Hydrological Institute ([SMHI](https://www.smhi.se/en/q/Stockholm/2673730)) and the UK Met Office Hadley Centre ([UKMO](https://www.metoffice.gov.uk/weather/climate/met-office-hadley-centre/index)) to develop and test a Multi-Envelope s-coordinate system ([Bruciaferri et al. 2018](https://link.springer.com/article/10.1007/s10236-018-1189-x)) in the Baltic sea.

This repository is still under construction. It has begun to be populated so that version tracking and branching can begin.

## Quick-start

```shell
git clone https://github.com/oceandie/BALTIC-MEs.git
cd BALTIC-MEs
conda env create -f pyogcm.yml
conda activate pyogcm
```

## Directories structure:

```
BALTIC-MEs
 `-- src/                            <-- source code
    |-- f90/                         <-- Fortran 90 code (mainly NEMO)
    |   `-- NEMO_4.0-HEAD_test_MEs/  <-- NEMO_4.0-HEAD_test_MEs UKMO branch rev 15806
    |       |-- arch/
    |       |-- cfgs/
    |       |-- ext/
    |       |-- mk/
    |       |-- src/
    |       |-- tests/
    |       `-- tools/               <-- tools_r4.0-HEAD_dev_MEs UKMO branch rev 15807
    `-- python/                      <-- Python3 code
        |-- envelopes/               <-- code to generate envelopes
        |   `-- templates_inp_files/ <-- template inlut files to generate AMM15, zco or sco envelopes
        `-- loc_area/                <-- code to generate localisation area
```
