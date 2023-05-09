# Historic-DICE

To results of "The Social Cost of Fossil and Industrial CO2 Emissions from 1950 to 2018" are build on 3 steps:
1) Calibration, 2) SCC Calculations, 3) IW Calculations. 

Each folder contains a separate readme file to give further instructions. However, none of the folders requires that the previous folder has been started, e.g., the IW calculations can be carried out without the Calibration or SCC_Calculations.

01_Calibration: The folder contains all input data,  the files for the historical spin-up and the calibration of the exogenous parameters and the other components of integrated assessment framework. As part of the spin-up, also the calculation of the share of atmospheric CO2 in the atmosphere, attributed to the CO2 emissions of an individual country is carried out. 

02_SCC_Calculations: Includes all AMPL files, distinguished for the calibrated and expert-based SDR. In each folder is a model file (which is identical) and then different .run files for the different impact functions.

03_IW_calculations: Includes all Matlab files to calculate the global IW assessment and the country IW assessment, both for gross and net investments and the calibrated and expert-based SDR.
