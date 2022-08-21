# Generating envelopes
`cd envelopes/`

1. #### Standard z-coordinates
`python generate_envelopes.py templates_inp_files/MEs_1env_zco.inp`

2. #### Standard s-coordinates
`python generate_envelopes.py templates_inp_files/MEs_1env_sco.inp`

3. #### UKMO [AMM15](https://github.com/JMMP-Group/CO_AMM15) Multi-Envelope model - NO local optimisation
`python generate_envelopes.py templates_inp_files/MetOffice_AMM15/amm15_MEs_2env_r01-003.inp`

4. #### UKMO [AMM15](https://github.com/JMMP-Group/CO_AMM15) Multi-Envelope model - 2 iterations of local optimisation
`python generate_envelopes.py templates_inp_files/MetOffice_AMM15/amm15_MEs_2env_r01-003_v2.inp`
