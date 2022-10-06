Code and data related to the paper 'Assessing the effectiveness of environmental sampling for surveillance of foot-and-mouth disease virus in a cattle herd'.

MATLAB REQUIREMENTS AND CODE:
The scripts/functions were run using Matlab version 2020b and require the Statistics and Machine Learning and Parallel Computing toolboxes. 
However, they can be easily adapted to run without the Parallel Computing toolbox by changing the "parfor" loop in Detection_contours.m to a "for" loop.

RandomPosterior_Model.m - Calls the Transmission_model function multiple times using random posterior distributions (from the five 2007 IPs).
Transmission_model.m - Runs the transmission model with both direct and environmental transmission.
Direct_model.m - Runs the transmission model with direct transmission only.
Env_model.m -Runs the transmission model with environmental transmission only.

Load_parameters.m - Loads either the prior or posterior distributions for each parameter.
Viral_profile.m - Computes the viral profile for each animal.

Envdynamics.m - Computes the probability of a sample detecting FMDV in different areas of the environment and generates Fig 1.

After running RandomPosterior_Model.m, run one of the following to generate figures:
Detection_Contours.m - Generate contour plots for detection time and thetas (e.g. Figs 2 and 3).
FreedomFromInfection.m - Generate results relating to confidence of freedom of infection (e.g. Figs 4-7 and Tables 2 and 3).

Data folder:
EnvSampProb - Posterior data for the probability of detection.
EnvSampProb_PCR - Posterior data for the probability of detection from PCR.
EnvTrans_SimpleModelToo4_MCMCSamples.mat - Prior data for environmental contamination and transmission parameters.
EnvironmentalTranmsissionData.mat - Prior data for environmental contamination and transmission parameters.
Posteriordistributions.mat - Posterior data for all five farms from the ABC-SMC algorithm.
