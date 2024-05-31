library(dplyr)

#=========================
#Dealing with Fc standards
#=========================

# Get the subfolder path
subfolder <- file.path(getwd(), 
                       "MBPdG10")

# Create a list for storing Fc half wave potentials (E1/2)
Fc_half_wave_potentials <- list()

# Get file name matching the pattern "*_Fc.csv" in the MBPdG10 folder
Fc_file_list <- list.files(subfolder, 
                           pattern = "*_Fc.csv", 
                           full.names = TRUE)

# Create an empty dataframe to store the results
Fc_E_half <- data.frame(FilePattern = character(),
                        HalfWaveAverage = numeric(),
                        stringsAsFactors = FALSE)

# Loop through the file list and perform calculations on each file
for (file_name in Fc_file_list) {
  Fc_temp_df <- read.csv(file_name,
                         fileEncoding="latin1") %>% # Read data into temporary dataframe
    rename(Potential = "WE.1..Potential..V.", 
           Current = "WE.1..Current..A.") %>% # Rename the columns for working electrode
    filter(Scan == 2,
           between(Potential, -0.1, 0.7)) %>% # Filter to work only on scan 2 and in values Fc comes
    mutate(half_wave = (Potential[which.max(Current)] 
                        + Potential[which.min(Current)]) / 2) # Find half wave Fc
  
  # Calculate the average of half wave values. Needed because it repeats
  # entries w/mutate. There is only one value so the avg is the same as that
  # one value.
  half_wave_avg <- mean(Fc_temp_df$half_wave, 
                        na.rm = TRUE)
  
  # Get the ligand from the file name
  file_pattern <- gsub("_Fc.csv", 
                       "", 
                       basename(file_name))
  
  # Add the ligand name and half wave average to the dataframe
  Fc_E_half <- rbind(Fc_E_half, 
                     data.frame(Ligand = file_pattern,
                                Fc_E = half_wave_avg,
                                stringsAsFactors = FALSE))
} # End of for loop

# The values seen for Fc as an internal standard can vry quite a lot with some
# ligands, therefore these must be excluded, To do so, we are applying an 
# arbitrary cutoff of 10% difference to the external ferrocene standard
ext_value <- Fc_E_half$Fc_E[Fc_E_half$Ligand == "Ext"] # This finds the value of the external Fc ref in the dataframe
Fc_E_filtered <- Fc_E_half[abs(Fc_E_half$Fc_E - ext_value) <= 0.1 * ext_value, ] # This filters out values which are more than 10% different

# Calculate the average of the remaining values, so you are combining both 
# the external and internal standard Fcs here
Fc_E_avg <- mean(Fc_E_filtered$Fc_E, 
                 na.rm = TRUE)