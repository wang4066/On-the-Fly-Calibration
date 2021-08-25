# On-the-Fly-Calibration


The four files that contain important functions are:

mpwi. R: Maximum posterior weighted information criterion, which uses progressive restrictive MPWI adaptive item assignment

par_update.R: Bayesian update rule to dynamically update both item and person parameters in the 2PNO model

oc1.R: Traditional online calibration method: multiple EM cycles (MEM)

oc2.R: Traditional online calibration method: multiple EM cycles (MEM), the difference is in line 65-67, with the addition of
 if(M==0){
        liq[i,] <- NA
      }else{
 
This is used to handle a special case in which a learner receives no operational items, which will happen when the session length 
is 20 and there are only 20% operational items in the bank. 

Because MEM runs much slower than par_update.R (the Bayesian moment matching update rule), we used oc1.R whenever possible and oc2.R only
in cases when there are learners who received no operational items.

To replicate the four simulation conditions, the following four files can be run in their entirety. 

Simulation_mpwi.R: Which corresponds to the MPWI results in the paper 

Simulation_md.R: Which corresponds to the match difficulty results in the paper

Simulation_rd.R: Which corresponds to the random selection results in the paper

Simulation_OC.R: Which corresponds to the online calibration results in the paper


