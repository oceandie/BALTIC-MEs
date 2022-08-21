## Overview
This is UKMO code to generate a domain configuration file (aka `domain_cfg.nc`) for models employing a Multi-Envelope s-coordinate system ([Bruciaferri et al. 2018](https://link.springer.com/article/10.1007/s10236-018-1189-x)) and based on NEMO v4.x .

The code includes the following UKMO branches:

1. **NEMO_4.0-HEAD_test_MEs/** is a copy of [NEMO_4.0-HEAD_test_MEs](https://code.metoffice.gov.uk/svn/nemo/NEMO/branches/UKMO/NEMO_4.0-HEAD_test_MEs) @ 15806

2. **NEMO_4.0-HEAD_test_MEs/tools** is a copy of [tools_r4.0-HEAD_dev_MEs](https://code.metoffice.gov.uk/svn/nemo/NEMO/branches/UKMO/tools_r4.0-HEAD_dev_MEs) @ 15807

---

## Compiling the DOMAINcfg tool
The domain configuration file is generated via the NEMOv4 DOMAINcfg tool. In the UKMO supercomputer (Monsoon) the F90 DOMAINcfg tool is compiled as:

1. `cd NEMO_4.0-HEAD_test_MEs/tools`

2. `./maketools -m XC40_METO_debug -n DOMAINcfg` where _XC40_METO_debug_ refers to `NEMO_4.0-HEAD_test_MEs/arch/UKMO/arch-XC40_METO_debug.fcm`

In order to compile the DOMAINcfg tool with a machine which is not listed in `NEMO_4.0-HEAD_test_MEs/arch/` an adhoc arch file configuration spcific for the new machine should be created.  
