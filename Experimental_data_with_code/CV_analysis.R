library(dplyr)

# Load the data that is in the file which extract the 
# value for E1/2 Fc to nomralise data
source("Fc_standard.R")
source("SW_analysis.R")

#===================
#===================
#===================
# Cyclic Voltammetry
#===================
#===================
#===================

#===================
#===================
# Scan rate studies
#===================
#===================

#=====================================
# File selection for scan rate studies
#=====================================
# Find the files with the CV data 
CV_list <- list.files(file.path("./MBPdG10"),
                           pattern = ".*_isol.*\\.csv$",
                           full.names = TRUE)

# Filter for the first peak isolated and scanrate
CV_isol_list <- CV_list[!grepl("_all_",
                               CV_list)] 

# Filter for the all peaks and scanrate
CV_all_list <- CV_list[grepl("_all_",
                             CV_list)]

# Filter such that a dataset is kept for each ligand. If a file set exists for 
# isolating the first oxidation peak at different scan rates, that is kept 
# preferrentially over that containing all oxidation peaks
filter_unique_ligands_scanrate <- function(file_list) {
  # Extract unique ligands
  unique_ligands <- unique(gsub("(.*)_isol.*\\.csv$", 
                                "\\1", 
                                file_list))
  
  # Initialize a vector to store the selected files
  selected_files <- character()
  
  # Loop through each unique ligand
  for (ligand in unique_ligands) {
    # Check if both "isol" and "isol_all" versions exist for the ligand
    if (any(grepl(paste0("^", 
                         ligand, 
                         "_isol_"), 
                  file_list)) && 
        any(grepl(paste0("^", 
                         ligand, 
                         "_isol_all_"), 
                  file_list))) {
      # Priority to "isol" version
      selected_files <- c(selected_files, 
                          file_list[grepl(paste0("^", 
                                                 ligand, 
                                                 "_isol_"), 
                                          file_list)])
    } else {
      ## If only one version exists or neither, add whatever version is present
      selected_files <- c(selected_files, 
                          file_list[grepl(paste0("^", 
                                                 ligand, 
                                                 "_isol"), 
                                          file_list)])
      
    }
  }
  
  return(selected_files)
}

CV_scanrate_ligands <- filter_unique_ligands_scanrate(CV_list)

#====================================
# Data cleaning for scan rate studies
#====================================

# Create empty dataframes for CVs of first isolated oxidative peak w/scan rate
CV_blank_df <- data.frame(Potential = numeric(),
                           Current = numeric(),
                           Scan_rate = numeric(),
                           Ligand = character())

# Function for extracting ligand name & scan rate from file name and setting the
# potential to be referenced against Fc. Accepts a list of file names and 
# a dataframe into which the data should be loaded as arguments
analyse_CV_scanrate <- function(file_list,df)
{
  for (file_name in file_list) {
  # Get the ligand from the file name
  L_name <- gsub("(.*)_isol.*\\.csv$", 
                 "\\1",
                 basename(file_name))
  
  Scanrate <- gsub("[^0-9.]",
                    "",
                   gsub(".*_isol_(.*)\\.csv$", 
                 "\\1",
                 basename(file_name))
                )
  
  # Create a temporary dataframe to load in and manipulate data
  CV_temp_df <- read.csv(file_name,
                         fileEncoding="latin1") # Read data into temporary dataframe
  
  # Check if column names need adjustment
  if (!all(c("Potential", "Current") %in% colnames(CV_temp_df))) {
   # Adjust column names
    colnames(CV_temp_df) <- c("Potential", 
                              "Uncorrected Current", 
                              "Current", 
                              colnames(CV_temp_df)[-(1:3)])
  }
  
  # Filter and rearrange the columns
  CV_temp_df <- CV_temp_df[, 
                           c("Potential", 
                               "Current",
                                "Scan")]
  
  # Adjust potential and add ligand name column
  CV_temp_df <- CV_temp_df %>%
    mutate(Potential = (Potential - Fc_E_avg)) %>% # Adjust potential to be vs Fc E1/2
    mutate(Ligand = L_name) %>% # Add a column which gives the ligand name
    mutate(Scan_rate = as.numeric(Scanrate)) # Add a column which gives scanrate
  
  # Append the temporary modified data to a single dataframe with all ligand data
  df <- rbind(df,
              CV_temp_df)
  } # End of for loop
  return(df)
} # End of function

# Create a dataframe which contains the data only for isolating the first peak
# at different scan rates
CVs_isolated_peaks_scanrate_df <- analyse_CV_scanrate(CV_isol_list,
                                    CV_blank_df)

# Create a dataframe which contains the data only for all oxidation peaks at
# different scan rates
CVs_all_peaks_scanrate_df <- analyse_CV_scanrate(CV_all_list,
                               CV_blank_df)

# Create a dataframe which ensures that scanrate dependence data is included for
# every ligand measured, irrespective of whether it corresponds to a scan rate
# study of an isolated peak or all oxidative peaks
CVs_all_ligands_scanrate_df <- analyse_CV_scanrate(CV_scanrate_ligands,
                                                 CV_blank_df)

#==========================================
# Peak picking in CVs for scan rate studies
#==========================================

# This function finds CV peaks for different scan rates. It uses the SW derived
# first oxidation potentials as a reference to restrict the window of searching
# for maximum current to +-100 mV from that value.
# Arguments: df1 - the dataframe containing the CV data at different scanrates
# df2 - the dataframe containing the reference SW data
CV_peak_find_scanrate <- function(df1, df2) {
  result_df <- data.frame()
  
  # Iterate over unique Ligand values in df1
  for (ligand in unique(df1$Ligand)) {
    
    # Iterate over unique Scan_rate values for each Ligand
    for (scan_rate in unique(df1$Scan_rate[df1$Ligand == ligand])) {
      
      # Extract the corresponding potential value from df2
      potential <- df2$Potential[df2$Ligand == ligand]
      
      # Check if there is a matching Ligand in df2
      if (length(potential) > 0) {
        # Calculate the search range
        search_range <- c(potential - 0.1, 
                          potential + 0.1)
        
        # Filter df1 based on Ligand, Scan_rate, and search range
        filtered_df <- df1 %>%
          filter(Ligand == ligand, 
                 Scan_rate == scan_rate,
                 Potential >= search_range[1],
                 Potential <= search_range[2],
                 Scan == 2)
        
        # Find the row with the maximum current within the search range
        max_current <- filtered_df[which.max(filtered_df$Current), 
                                   c("Ligand", 
                                     "Potential",
                                     "Current",
                                     "Scan_rate",
                                     "Scan")]
        
        # Add the result to the output dataframe
        result_df <- bind_rows(result_df, 
                               max_current)
      }
    }
  }
  
  return(result_df)
}

# Use function to get peaks at different scanrates for all ligands; arrange
# by ligand and in increasing order of scanrate for sanity check of values.
ligands_scanrate <- CV_peak_find_scanrate(CVs_all_ligands_scanrate_df,
                                first_oxidation_potential_all) %>% # Use above function to get CV peaks
  arrange(Ligand,
          Scan_rate) %>% # Arrange it by ligand and then in ascending order of scan rate
  mutate(Scan_rate_log = log10(Scan_rate/1000)) %>% # Create a column with log(scan rate)
  mutate(Scan_rate_sqrt = sqrt(Scan_rate))

# Get the coefficients of the models in the above dataframe
coefficients_Potential_vs_logScanrate <- ligands_scanrate %>%
  group_by(Ligand) %>% # Ensures you apply following operations to each ligand
  summarise(model_log = list(
                lm(
                  Potential ~ Scan_rate_log))) # Linear model between potential and log(scan rate) in V/s for each ligand

# Coefficients for the linear models between current and square root of scan rate
coefficients_Current_vs_sqrtScanrate <- ligands_scanrate %>%
  group_by(Ligand) %>%
  summarise(model_sqrt = list(
    lm(
      Current ~ Scan_rate_sqrt)))

#=======================
#=======================
# Scan direction studies
#=======================
#=======================

#==========================================
# File selection for scan direction studies
#==========================================

# Find the files with the CV data for the first peak isolated
CV_full_list <- list.files(file.path("./MBPdG10"),
                           pattern = ".*_full.*\\.csv$",
                           full.names = TRUE)

#=========================================
# Data cleaning for scan direction studies
#=========================================
# Create empty dataframe for CVs of first isolated 
# oxidative peak w/scan rate
CV_df_blank <- data.frame(Potential = numeric(),
                          Current = numeric(),
                          Scan_direction = character(),
                          Ligand = character())

# Function for extracting ligand name & scan direction from file name and setting the
# potential to be referenced against ferrocene. Accepts a list of file names and 
# a dataframe into which the data should be loaded as arguments
analyse_CV_scandir <- function(file_list,df)
{
  for (file_name in file_list) {
    # Get the ligand from the file name
    L_name <- gsub("(.*)_full.*\\.csv$", 
                   "\\1",
                   basename(file_name))
    
    Scandir <- gsub(".*_full_(.*)\\.csv$", 
                          "\\1",
                          basename(file_name))
    
    # Create a temporary dataframe to load in and manipulate data
    CV_temp_df <- read.csv(file_name,
                           fileEncoding="latin1") # Read data into temporary dataframe
    
    # Check if column names need adjustment
    if (!all(c("Potential", "Current") %in% colnames(CV_temp_df))) {
      # Adjust column names
      colnames(CV_temp_df) <- c("Potential", 
                                "Uncorrected Current", 
                                "Current", 
                                colnames(CV_temp_df)[-(1:3)])
    }
    
    # Filter and rearrange the columns
    CV_temp_df <- CV_temp_df[, 
                             c("Potential", 
                               "Current",
                               "Scan")]
    
    # Adjust potential and add ligand name column
    CV_temp_df <- CV_temp_df %>%
      mutate(Potential = (Potential - Fc_E_avg)) %>% # Adjust potential to be vs Fc E1/2
      mutate(Ligand = L_name) %>% # Add a column which gives the ligand name
      mutate(Scan_direction = Scandir) # Add a column which gives scanrate
    
    # Append the temporary modified data to a single dataframe with all ligand data
    df <- rbind(df,
                CV_temp_df)
  } # End of for loop
  return(df)
} # End of function

# Apply above function to existing CV data
CV_dir <- analyse_CV_scandir(CV_full_list,
                             CV_df_blank)

# This function finds CV peaks without scan rates. It provided
# estimated potentials as a reference to restrict the window of searching
# for maximum current to +-100 mV from that value.
# Arguments: df1 - the dataframe containing the CV data at different scanrates
# df2 - the dataframe containing the reference SW data
CV_peak_find <- function(df1, df2) {
  result_df <- data.frame()
  
  # Iterate over unique Ligand values in df1
  for (ligand in unique(df1$Ligand)) {
      
      # Extract the corresponding potential value from df2
      potential <- df2$Potential[df2$Ligand == ligand]
      
      # Check if there is a matching Ligand in df2
      if (length(potential) > 0) {
        # Calculate the search range
        search_range <- c(potential - 0.1, 
                          potential + 0.1)
        
        # Filter df1 based on Ligand, Scan_rate, and search range
        filtered_df <- df1 %>%
          filter(Ligand == ligand,
                 Potential >= search_range[1],
                 Potential <= search_range[2])
        
        # Find the row with the maximum current within the search range
        max_current <- filtered_df[which.max(filtered_df$Current), 
                                   c("Ligand", 
                                     "Potential")]
        
        # Add the result to the output dataframe
        result_df <- bind_rows(result_df, 
                               max_current)
      }
    }
  
  return(result_df)
}

# Find reduction peaks for ligands which should have genuine ones, based on scan
# direction studies
# Create a dataframe with estimated values for this reduction
approx_E_red <- data.frame("Potential" = c(-1.8,
                                           -1.7,
                                           -1.8),
                           "Ligand" = c("PtBu3",
                                        "PtBu2Np",
                                        "PtBu2Cy"))

# Filter the scan direction data to include only the ligands which have a 
# reductive peak in the first negative scan and filter for the second scan
CV_dir_reductive <- CV_dir %>%
  filter(Ligand %in% c("PtBu3","PtBu2Np","PtBu2Cy"),
         Scan_direction == "pos",
         Scan == 1)

# Use the function to find the peaks
reductive_peaks_CV <- CV_peak_find(CV_dir_reductive,
                                   approx_E_red) %>%
  rename("Reduction Potential" = "Potential")
