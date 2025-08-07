# RE

Data and code for the publication:

Metal Electronics Mediate Steric and Coordination Number Effects on Palladium (II) C–X Reductive Elimination

by

Michael K. Bogdos, Sven Roediger, Florian Ruepp, Nathalie A. V. Rowlinson, Patrick Müller, Fabio Masero, Jan Hübscher, Bill Morandi

DFT directory contains:
  Energy_decomposition_analysis: 
  - DFT output files for energy decomposition analysis. Nomenclature: L_state.out; L - ancillary ligand on F4-NClMs-OMe-L palladacycle. state - ground state (gs) or transition state (ts).
  
  Lit_complexes: 
  - xyz files for the complexes and fragments used for parametrization of the literature dataset.

  Ni_complexes_monodentate_phosphines:
  - xyz, output and Hessian files for the LNi(CO)3 complexes used for calculating Tolman parameters.
  
  Palladacycles: 
  - xyz and log files for ground states and transition states of the palladcycle reductive elimination. Nomenclature: (state)_(L).out; (L) - ancillary ligand on F4-NClMs-OMe-L palladacycle. (state) - ground state (gs) or transition state (ts).

  Rh_complexes_bidentate_phosphines: 
  - xyz, output and Hessian files for the LNi(CO)3 complexes used for calculating CO stretching frequency of Rh(H)(CO)L complexes.

experimental_data_without_code directory contains:
  NMR and IR raw data
  HRMS core facility reports

experimental_data_with_code directory contains:
Electrochemical data
  - raw data in MBPdG10
  - code used for data processing 
    (CV_analysis.R, SW_analysis.R, FC_standard.R)
  - code used for plotting
    (CV_plots.R, SW_plots.R, Eox_descriptors.R)

Kinetic data for palladacycle C-N reductive elimination
  - raw data for different backbones in subdirectory half_life_data
  - raw data for different ancillary ligands on F4-NClMs-OMe backbone in subdirectory eyring
  - data analysis and plotting for different backbones
    (backbone_half_lives.R)
  - data analysis and plotting for eyring analysis 
    (eyring.R)
  - error analysis for the various kinetics 
    (kinetics_errors.R)

Modeling
  - power analysis and seeds for reproducible simulations
    (generate_simulated_data.R, analyse_simulations.R, simulations_plots.R, seeds.txt)
  - simulations for the potential for endogenous error that scales with logk
    (relative_error_simulations.R)
  - model selection on palladacycle data for varying ancillary ligands
    (ligand_effect_analysis.R)
  - mediation analysis and path coefficients
    (mediation_analysis.R)

DFT data
  - VIP, AIP, Vbur etc. in subdirectory complex_descriptors
  - extracted benchmarked barrier data for palladacycles in subdirectory barriers

Literature data
  - data from the literature and DFT calculated values in subdirectory literature_data
  - energy decomposition analysis extracted data in subdirectory EDA

Package information for the code
  - lock file renv.lock contains all package names and versions
