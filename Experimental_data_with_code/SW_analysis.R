library(dplyr)
library(stringr)

# Load the data that is in the file which extract 
# the value for E1/2 Fc to nomralise data 
source("Fc_standard.R")

#========================
#========================
# Square Wave Voltammetry
#========================
#========================

#==========
# Functions
#==========
# Find the files with the square wave data
SW_file_list <- list.files(file.path("./MBPdG10"), 
                           pattern = "SW_.*\\.csv", 
                           full.names = TRUE)

# A fucntion which loops through a file 
# list and combines & manipulates all the SW data
# the argument just requires a list 
# of files with the file name structured as 
# SW_LIGAND.csv to function
SW_comb <- function(file_list){
  
  # Create empty dataframe which will 
  # be filled with the rest
  SWs_all <- data.frame()
  
  for (file_name in file_list) {
  
  # Get the ligand from the file name
  L_name <- gsub("SW_(.*)\\.csv", 
                 "\\1",
                 basename(file_name))
  
  # Create a temporary dataframe to 
  # load in and manipulate data
  SW_temp_df <- read.csv(file_name,
                         check.names = F) 
  
  # Check if column names need adjustment
  if (!all(c("Potential", 
             "Current") %in% colnames(SW_temp_df))) {
    
   # Adjust column names
    colnames(SW_temp_df) <- c("Potential", 
                              "Uncorrected Current", 
                              "Current", 
                              colnames(SW_temp_df)[-(1:3)])
  } # End of if clause
  
  # Filter and rearrange the columns
  SW_temp_df <- SW_temp_df[, 
                           c("Potential", 
                               "Current")]
  
  # Adjust potential and add ligand name column
  SW_temp_df <- SW_temp_df %>%
    mutate(Potential = (Potential - Fc_E_avg)) %>%
    mutate(Ligand = L_name)
  
  # Append the temporary modified data to a 
  # single dataframe with all ligand data
  SWs_all <- rbind(SWs_all,
                   SW_temp_df) 
  
  } # End of for loop
  
  # Final function output
  return(SWs_all)
}

# Make a function which picks peaks out of a 
# combined dataframe of Square Wave data
# The arguments that the function needs are: 
# df = dataframe (usually outputted by SW_comb)
# V = voltage window in which you want to look for peaks 
# I = current threshold for looking for peaks,
# no_peaks = the number of peaks you want to detect, 
# 1 will give you the first peak found
# starting from negative potentials and 
# going to more positive. The variable I should not be set 
# too low as the baseline is small 
# but technically sinusoidal so it will 
# give peaks detected where
# there are none really.
# Tldr; SW_peaks(data,
# Voltage_threshold,
# Current_threshold,number of peaks to be detected)
SW_peaks <- function(df,V,I,no_peaks){
  
  # Append data to peak potential data frame
  peak_potential_df <- df %>%
    filter(Potential >= V & Current >= I) %>%
    filter(row_number() %in% (which(
      diff(sign(diff(Current))) < 0) + 1)) %>%
    arrange(Potential) %>%
    group_by(Ligand) %>%
    slice_head(n = no_peaks) %>%
    ungroup() %>%
    arrange(Potential)
  
  # Final output of the function
  return(peak_potential_df)
  
}

#=========
# Datasets
#=========
# Square wave data
# All palladacycles - add a column that classifies 
# the coord no
all_SW <- SW_comb(SW_file_list) %>%
  filter(!(str_detect(Ligand,
                    "_free"))) %>%
  mutate(L_type = ifelse(
    Ligand %in% c("BIPHEP", 
                "dcype", 
                "dppe", 
                "Fdppe",
                "dppm",
                "dmpe",
                "dppbz"),
    "bidentate",
    "monodentate"))

# All free ligands
all_free_L_SW <- SW_comb(SW_file_list) %>%
  filter(str_detect(Ligand,
                    "_free"))

# Extracted peaks data
# Get the first oxidation potential 
# All the ligands
first_oxidation_potential_all <- SW_peaks(all_SW,
                                          0.4,
                                          5e-7,
                                          1)

# Bidentate
first_oxidation_potential_bidentate <- SW_peaks(all_SW %>%
                                                  filter(
                                        L_type == "bidentate"),
                                                0.5,
                                                1e-6,
                                                1)

# Monodentate
first_oxidation_potential_monodentate <- SW_peaks(all_SW %>%
                                                    filter(
                                      L_type == "monodentate"),
                                                  0.5,
                                                  1e-6,
                                                  1)

# All free ligands
first_ox_free_L <- SW_peaks(all_free_L_SW,
                            0.3,
                            1e-6,
                            1)

# Get the first 4 peaks for all the data - this 
# should find them, there should be 2 or 3
# All ligands
first_three_oxidation_potentials_all <- SW_peaks(all_SW,
                         0.5,
                         1e-6,
                         4)

# Bidentate
first_three_oxidation_potentials_bidentate <- SW_peaks(
                                                all_SW %>%
                                                  filter(
                                      L_type == "bidentate"),
                                                0.5,
                                                1e-6,
                                                4)

# Monodentate
first_three_oxidation_potentials_monodentate <- SW_peaks(
                                                  all_SW %>%
                                                    filter(
                                    L_type == "monodentate"),
                                                  0.5,
                                                  1e-6,
                                                  4)
