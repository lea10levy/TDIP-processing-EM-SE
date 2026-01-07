# TDIP-processing-EM-SE
MATLAB scripts for processing TDIP data including early time EM noise removal and subsequent decay fitting with stretched exponential

By Léa Lévy 
2026-01-07

Starting files = mat files for logging TDIP data and text files (.tx2 extension) for surface and cross-borehole TDIP data. All data are found in the “input data” folder.

1) PROCESSING
Krafla_EM_SEfit_v3.m works with surface TDIP data. It calls .tx2 files in Data repository\Input data\Krafla surface TDIP data\
Hvede_EM_SEfit_v3.m works with cross-borehole TDIP data. It calls .tx2 files in Data repository\Input data\Hvedemarken cross borehole TDIP data\
Nesjavellir_EM_SEfit_v3.m works with wireline logging data from the QL40-IP tool from ALT. It calls .mat files in Data repository\Input data\Nesjavellir logging TDIP data\

The three scripts are relatively similar. They load data and carry out the seven steps of the processing procedure described in the paper. 
They plot the results (as shown in Fig. 3 of the paper) and save the results as a .mat file which can be retrieved in 'Processing output files' folder.

2) INVERSION AND PLOTTING
tx2_to_tx2_processed_Krafla_Hvede_Nes.m converts the processed dataset from a .mat file to an inversion-ready .tx2 file (the subsequent converison to the .qui format is carried out through EEMstudio)
plot_qui_quo_decays.m plots the TDIP decays, both data and forward, to evaluate the fitting qsuality after inversion
