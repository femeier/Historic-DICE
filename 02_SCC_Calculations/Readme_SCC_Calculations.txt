For the calculation of the SCC timeseries, the IAM described in Section 3 of SI_B has been implemented in AMPL. 
The model file is DICE2016_Fair_Geoffroy_SCC.mod.

The file is called by a .run file whereby we have three different run files, one for each impact function each. 
R_DICE_SCC_Nordhaus.run
R_DICE_SCC_Sterner.run
R_DICE_SCC_Weitzman.run

The model is solved by using the Knitro Solver.

There are two folders, one for the calibrated SDR and one for the expert-based SDR.

Each folder contains a .run file, called "master_call.run". Starting this routine will start in the corresponding folder the three input files and all results will be calculated and written as .csv files to the folder. 
Note that the folders already contain the result files, i.e, starting the master_call.run file or an individual run file will overwrite the result files. 