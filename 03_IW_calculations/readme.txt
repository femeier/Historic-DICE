The folder contains three matlab files which perform the IW calculations on the global level (CWB_global_calibrated.m and CWB_global_expert.m) and on the country level (Country_results.m).

The folder contains also a subfolder, "Data" which has the background data. This is for information only. All input files are in "Input_data_raw_sorted.xlxs". 

The files for the global IW calculations have to be excuded first since they provide input for the country calculations. 

However, the current version already contains all output files (and hence inputfiles for subsequent calculations) which means that alo the country_results.m file can be started without running the global files first. 