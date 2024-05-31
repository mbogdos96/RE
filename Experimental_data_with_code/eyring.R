library(dplyr)
library(ggplot2)
library(nls2)
library(purrr)
library(patchwork)

source("SW_plots.R")
source("common_functions.R")

#  Read in csvs for subfolder
# Get a list of all CSV files in the eyring folder
csv_files <- list.files(path = "./eyring", 
                        pattern = "\\.csv$", 
                        full.names = TRUE)

import_eyring <- function(files){
  # Initiate list
  df_list <- list()
  
  for (file in files) {
    # Extract the file name without extension
    file_name <- tools::file_path_sans_ext(basename(file))
    
    # Read the CSV file into a data frame
    df <- read.csv(file)
    
    # Rename columns
    colnames(df) <- gsub("time\\.\\.\\.s", "time_s", 
                         gsub("c\\.Pd\\.Complex\\.\\.\\.\\.M", "Pd_molar", 
                              colnames(df)))
    
    # Covert time to numerics
    df <- df %>%
      mutate(time_s = as.numeric(time_s))
    
    # Store the data frame in the list
    df_list[[file_name]] <- df
  }
  
  return(df_list)
}

eyring_data_list <- import_eyring(csv_files)

# Function to fit first order rate constant to all data in a list 
# of dataframes and create ggplot object for each fit
fit_k_eyring <- function(df_list) {
  # Initialize an empty dataframe for output
  df_out <- data.frame(Ligand = character(), 
                       Backbone = character(), 
                       Temperature_K = numeric(), 
                       rate_constant = numeric())
  
  # Initialize an empty dataframe for predicted values
  df_predicted <- data.frame(Ligand = character(), 
                             Backbone = character(), 
                             time_s = numeric(),
                             Pd_predicted = numeric())
  
  # Loop through each dataframe in the list
  for (df in seq_along(df_list)) {
    # Fit the model
    fit <- nls(Pd_molar ~ A * exp(-B * time_s) + C, 
               data = df_list[[df]], 
               start = list(A = 1e-2,
                            B = 1e-4,
                            C = 1e-3),
               lower = c(0,0,0),
               algorithm = "port")
    
    # Extract values from df and predict Pd_molar 
    # values based on range of values
    # of time in original df
    sim_time_s <- seq(range(df_list[[df]]$time_s)[1],
                      range(df_list[[df]]$time_s)[2],
                      by = 10)
    predicted <- predict(fit, data.frame(time_s = sim_time_s))
    
    # Create dataframe for predicted values
    sim_df <- data.frame(Ligand = rep(df_list[[df]]$Ligand[1], 
                                      length(sim_time_s)),
                         Backbone = rep(df_list[[df]]$Backbone[1], 
                                        length(sim_time_s)),
                         time_s = sim_time_s,
                         Pd_predicted = predicted)
    
    # Merge the predicted df of the loop with the one outside the loop
    df_predicted <- bind_rows(df_predicted,
                              sim_df)
    
    # Extract ligand, backbone, T and coefficients
    df_row <- df_list[[df]] %>%
      slice(1) %>%
      dplyr::select(Ligand, 
                    Backbone, 
                    Temperature_K) %>%
      mutate(Pd_initial_molar = coef(fit)["A"],
             rate_constant = coef(fit)["B"],
             baseline_correction = coef(fit)["C"])
    
    # Append L,backbone, T, coef to output dataframe
    df_out <- bind_rows(df_out, 
                        df_row)
  }
  
  # Return output dataframe and list of ggplot objects
  return(list(Rate_constants_eyring = df_out,
              Predicted_values_eyring = df_predicted))
}

# Apply function to get rate constants and plots
eyring_parameters_and_predictions <- fit_k_eyring(eyring_data_list)

# Separate out the df used for the table in the md and fix naming
eyring_summary <- eyring_parameters_and_predictions$Rate_constants_eyring %>%
  mutate(rate_constant = unname(rate_constant),
         Pd_initial_molar = unname(Pd_initial_molar),
         baseline_correction = unname(baseline_correction),
         Barrier_kcal_mol = -8.314*Temperature_K/4184
         *log(rate_constant/(2.1e10*Temperature_K)))

# Function for creating plots and fits
eyring_plots <- function(pred_df,
                         df_list){
  # Initialize an empty list for ggplot objects
  ggplot_list <- list()
  
  for (df in seq_along(df_list)){
    # Fix the labels in each df
    temp_df <- data.frame("Ligand" = as.character(
      df_list[[df]]$Ligand[1]),
                          "Backbone" = as.character(
                            df_list[[df]]$Backbone[1])) %>%
    mutate(Backbone_fixed = sapply(Backbone,
                                   replace_backbone_labels),
           Ligand_fixed = sapply(Ligand,
                                 replace_ligand_expression_og),
           Fixed_name = paste(Backbone_fixed,
                              "*-",
                              Ligand_fixed))
    
    # Create ggplot object
    plot <- ggplot() +
      labs(title = parse(text = temp_df$Fixed_name),
        x = "Time (s)", 
        y = "Concentration (M)") +
      geom_point(data = df_list[[df]],
                 aes(x = time_s,
                     y = Pd_molar),
                 size = 3) +
      geom_line(data = filter(pred_df,
                              Ligand == df_list[[df]]$Ligand[1],
                              Backbone == df_list[[df]]$Backbone[1]),
                aes(x = time_s,
                    y = Pd_predicted),
                linewidth = 1,
                colour = "red") +
      themething
    
    # Name the ggplot object
    plot_name <- paste0(df_list[[df]]$Ligand[1], 
                        "_", 
                        df_list[[df]]$Backbone[1])
    
    # Add ggplot object to list
    ggplot_list[[plot_name]] <- plot
    
  }
  return(ggplot_list)
}

eyring_plot_list <- eyring_plots(
  eyring_parameters_and_predictions$Predicted_values_eyring,
  eyring_data_list)

eyring_plots_rm_NMe2 <- eyring_plot_list[-3]

# Wrap all eyring plots into one
all_eyring_plots <- wrap_plots(eyring_plots_rm_NMe2) +
  plot_layout(ncol = 3,
              nrow = 2,
              axis_titles = "collect") &
  theme(axis.text = element_blank(),
        plot.title = element_text(size = 12),
        axis.title = element_text(size = 10))

# Table for DFT benchmarking
dft_table_df <- eyring_summary %>%
    filter(Ligand %in% c("dppe","PtBu3")) %>%
    dplyr::select(Backbone, 
                  Ligand,
                  Barrier_kcal_mol) %>%
    rename(Experimental = Barrier_kcal_mol) %>%
  mutate("M06" = c(24.5,17.3),
         "Mo06-D3(0)" = c(24.3,17.3),
         "PBE0-D3(BJ)" = c(28.0,21.8),
         "wB97x-D" = c(30.5,22.7))