The excel file contains the various exogenous parameters.

Note that the calibration of eta and delta can be found in 02_SCC_Calculations\SCC_results_2022_calibrated_SDR since the calculations requires already the calibrated IAM. The IAM is optimized in the baseline specification (i.e. Nordhaus Impacts) with respect to investment and the deviation from the "optimal" capital stock against the "observed capita stock" is used as indicator for the goodness of fit. 

The AMPL run file is "R_DICE_Nordhaus_calibration.run" and the obtained capital stock is written to the file "out_combinations_eta_delta.csv".