# On-the-Fly-Calibration


The four files that contain important functions are:

mpwi. R: Maximum posterior weighted information criterion, which uses progressive restrictive MPWI adaptive item assignment

par_update.R: Bayesian update rule to dynamically update both item and person parameters in the 2PNO model

oc1.R: Traditional online calibration method: multiple EM cycles (MEM)

To replicate the four simulation conditions, the following four files can be run in their entirety. 

Simulation_mpwi.R: Which corresponds to the MPWI results in the paper 

Simulation_md.R: Which corresponds to the match difficulty results in the paper

Simulation_rd.R: Which corresponds to the random selection results in the paper

Simulation_OC.R: Which corresponds to the online calibration results in the paper


