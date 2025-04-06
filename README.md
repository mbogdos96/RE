# RE

Data and code for the publication:

Metal Electronics Mediate Steric and Coordination Number Effects on Palladium (II) C–X Reductive Elimination

by

Michael K. Bogdos, Sven Roediger, Florian Ruepp, Patrick Müller, Nathalie A. V. Rowlinson, Fabio Masero, Bill Morandi

DFT directory contains:
  DFT input and output files

experimental_data_without_code directory contains:
  NMR, IR and HRMS data

experimental_data_with_code directory contains:

 Electrochemical data
  - raw data in MBPdG10
  - code used for data processing 
    (CV_analysis.R, SW_analysis.R, FC_standard.R)
  - code used for plotting
    (CV_plots.R, SW_plots.R, Eox_descriptors.R)
 Kinetic data for palladacycle C-N reductive elimination
  - raw data for different backbones in half_life_data
  - raw data for different ancillary ligands on F4-NClMs-OMe backbone in eyring
  - data analysis and plotting for different backbones
    (backbone_half_lives.R)
  - data analysis and plotting for eyring analysis 
    (eyring.R)
  - error analysis for the various kinetics 
    (kinetics_errors.R)
 Modeling
  - power analysis and seeds for reproducible simulations
    (generate_simulated_data.R, analyse_simulations.R, simulations_plots.R, seeds.txt)
  - model selection on palladacycle data for varying ancillary ligands
    (ligand_effect_analysis.R)
  - mediation analysis and path coefficients
    (mediation_analysis.R)
 DFT data
  - VIP, AIP, Vbur etc. in complex_descriptors
  - extracted benchmarked barrier data for palladacycles in
 Literature data
  - data from the literature and DFT calculated values in literature_data
  - energy decomposition analysis extracted data in EDA
