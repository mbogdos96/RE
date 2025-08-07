library(dplyr)
library(purrr)

source("generate_simulated_data.R")

#===========================
# Individual metric analysis
#===========================
# Function for finding which model 
# each metric predicts is best
analyse_simulations <- function(df) {
  
  grouped_df <- df[[3]] %>%
    group_by(N_Obs, 
             Ground_Truth_Model,
             SD_Noise,
             Sim_no) %>%
    mutate(AICc_min = ifelse(AICc == min(AICc), 
                             Model, 
                             NA),
           BIC_min = ifelse(BIC == min(BIC), 
                            Model,
                            NA),
           MPIW_min = ifelse(Mean_Prediction_Interval_Width == 
                               min(Mean_Prediction_Interval_Width), 
                             Model, 
                             NA),
           SE_E_ox_min = ifelse(SE_Estimate_E_ox == min(SE_Estimate_E_ox), 
                                Model, 
                                NA),
           SE_Intercept_min = ifelse(SE_Estimate_Intercept == 
                                       min(SE_Estimate_Intercept), 
                                     Model, 
                                     NA),
           adj_R2_max = ifelse(adj_R2 == max(adj_R2), 
                                     Model, 
                                     NA),
           ANOVA_p_val_sel = case_when(
             all(ANOVA_log10p > log10(0.05)) ~ "1",
             all(ANOVA_log10p < log10(0.05)) ~ "3",
             ANOVA_log10p[1] > log10(0.05) 
             & all(ANOVA_log10p[2:3] < log10(0.05)) ~ "3",
             TRUE ~ "2")) %>%
    dplyr::select(N_Obs,
                  Ground_Truth_Model,
                  SD_Noise,
                  Sim_no,
                  AICc_min,
                  BIC_min,
                  MPIW_min,
                  SE_E_ox_min,
                  SE_Intercept_min,
                  adj_R2_max,
                  ANOVA_p_val_sel) %>%
    summarize(across(everything(), ~ first(na.omit(.)))) %>%
    ungroup
  
  return(grouped_df)
}

# Run single metric selection function on simulated data
analysed_simulations <- analyse_simulations(sim_results)

# Check if metrics select the ground truth model
calculate_individual_metric_truth <- function(df){
  df_out <- df %>%
    mutate(AICc_correct = ifelse(AICc_min == Ground_Truth_Model, 
                  T, 
                  F),
           BIC_correct = ifelse(BIC_min == Ground_Truth_Model, 
                                 T, 
                                 F),
           MPIW_correct = ifelse(MPIW_min == Ground_Truth_Model, 
                                 T, 
                                 F),
           SE_E_ox_correct = ifelse(SE_E_ox_min == Ground_Truth_Model, 
                                 T, 
                                 F),
           SE_Intercept_correct = ifelse(SE_Intercept_min == Ground_Truth_Model, 
                                    T, 
                                    F),
           adj_R2_correct = ifelse(adj_R2_max == Ground_Truth_Model, 
                                         T, 
                                         F),
           ANOVA_correct = ifelse(ANOVA_p_val_sel == Ground_Truth_Model, 
                                   T, 
                                   F))
  
  return(df_out)
}

# Run the truth finding function on the simulated data
binary_truth_df <- calculate_individual_metric_truth(
  analysed_simulations)

# Function for calculating probabilities over simulations for given levels of
# noise and ground truth
calc_probabilities_individual_metrics <- function(df){

  # Find columns containing "_correct"
  correct_cols <- grep("_correct$", 
                       names(df), 
                       value = T)
  
  # Create a new dataframe for each unique value of "_correct" columns
  correct_df <- lapply(correct_cols, function(correct_col) {
    temp_df <- df %>%
      group_by(N_Obs,
               Ground_Truth_Model,
               SD_Noise) %>%
      summarise(quotient = sum(.data[[correct_col]] == T) / n())
    temp_df$correct_col <- correct_col
    return(temp_df)
  })
  
  # Combine the dataframes into one
  final_df <- do.call(rbind, correct_df)
  
  return(final_df)
}

# Run single metric probabilitity calculation function 
# on simulated data
prob_ind_metric <- calc_probabilities_individual_metrics(
  binary_truth_df)

#============================
# Interrogate ANOVA selection
#============================
# Function for finding which model each 
# ANOVA interrogation metric is best
analyse_ANOVA <- function(df) {
  
  grouped_df <- df[[3]] %>%
    group_by(N_Obs, 
             Ground_Truth_Model,
             SD_Noise,
             Sim_no) %>%
    mutate(ANOVA_std_sig = case_when(
      all(ANOVA_log10p > log10(0.05)) ~ "1",
      all(ANOVA_log10p < log10(0.05)) ~ "3",
      ANOVA_log10p[1] > log10(0.05) 
      & all(ANOVA_log10p[2:3] < log10(0.05)) ~ "3",
      TRUE ~ "2"),
           ANOVA_high_sig = case_when(
             all(ANOVA_log10p > log10(0.1)) ~ "1",
             all(ANOVA_log10p < log10(0.1)) ~ "3",
             ANOVA_log10p[1] > log10(0.1) 
             & all(ANOVA_log10p[2:3] < log10(0.1)) ~ "3",
             TRUE ~ "2"),
           ANOVA_low_sig = case_when(
             all(ANOVA_log10p > log10(0.001)) ~ "1",
             all(ANOVA_log10p < log10(0.001)) ~ "3",
             ANOVA_log10p[1] > log10(0.001) 
             & all(ANOVA_log10p[2:3] < log10(0.001)) ~ "3",
             TRUE ~ "2")) %>%
    dplyr::select(N_Obs,
                  Ground_Truth_Model,
                  SD_Noise,
                  Sim_no,
                  ANOVA_std_sig,
                  ANOVA_high_sig,
                  ANOVA_low_sig) %>%
    summarize(across(everything(), ~ first(na.omit(.)))) %>%
    ungroup
  
  return(grouped_df)
}

# Apply the ANOVA criteria selection to the 
# simulated data
ANOVA_analysis_df <- analyse_ANOVA(sim_results)

# Check if ANOVA conditions select ground truth
ANOVA_truth <- function(df){
  df_out <- df %>%
    mutate(std_sig_correct = ifelse(ANOVA_std_sig == Ground_Truth_Model, 
                                T, 
                                F),
           high_sig_correct = ifelse(ANOVA_high_sig == Ground_Truth_Model, 
                                 T, 
                                 F),
           low_sig_correct = ifelse(ANOVA_low_sig == Ground_Truth_Model, 
                                    T, 
                                    F))
  return(df_out)
}

# Apply ground truth check to ANOVA criteria
ANOVA_truth_df <- ANOVA_truth(ANOVA_analysis_df)

# Apply probability calculation function to ANOVA df
ANOVA_prob <- calc_probabilities_individual_metrics(
  ANOVA_truth_df)

#=================================
# Combinations of metrics analysis
#=================================
# Function for calculating probability 
# correct for combination of m - n metrics
probability_metric_combos <- function(df, 
                                      n_met, 
                                      n_sims){
  
  # Filter out simulation number column as 
  # it is not needed
  df <- df %>%
    dplyr::select(-Sim_no)
  
 # Specify the columns to group by
 columns_to_group <- c("N_Obs", 
  "Ground_Truth_Model", 
  "SD_Noise")

  # Find all other columns
 other_cols <- setdiff(names(df), 
                       columns_to_group)
 
 # df which contains rows where all metrics agree
 agree_df <- df %>%
   group_by(across(all_of(columns_to_group))) %>%
   rowwise() %>%
   filter(n_distinct(c_across(other_cols)) <= n_met) 
 
 # Function to calculate mode
 calculate_mode <- function(x) {
   ux <- unique(x)
   ux[which.max(tabulate(match(x, ux)))]
 }
 
 # df which contains rows where all metrics agree and are correct
 correct_df <- agree_df %>%
   filter(calculate_mode(c_across(everything())) 
          == Ground_Truth_Model) %>%
   ungroup() %>%
   dplyr::count(N_Obs,
                Ground_Truth_Model,
                SD_Noise) %>%
   dplyr::rename("n_corr" = "n")
 
 # count the columns which agree
 agreement <- agree_df %>%
   ungroup() %>%
   dplyr::count(N_Obs,
                Ground_Truth_Model,
                SD_Noise)
 
 # df with only combos of models and noise to populate with probs
 only_vals <- df  %>%
   dplyr::select(N_Obs,
                 SD_Noise,
                 Ground_Truth_Model) %>%
   unique()
 
 # combine dfs and calculate probs
 final_df <- list(only_vals,
                  agreement,
                  correct_df) %>%
   purrr::reduce(left_join,
                 by = c("N_Obs",
                        "Ground_Truth_Model",
                        "SD_Noise")) %>%
   mutate(prob_agree = ifelse(!is.na(n),
                              n/n_sims,
                              0),
          prob_correct = ifelse(n != 0 & !is.na(n_corr),
                                n_corr/n,
                                0),
          abs_prob_correct = ifelse(n != 0 & !is.na(n_corr),
                                n_corr/n_sims,
                                0))
          
 return(final_df)}

# Function for applying the probability function 
# to a list of  numbers which determines how many 
# are left out
leave_how_many_out <- function(df,numbers,n_sims){
  
  # Initialize df list
  final_df <- list()
  
  # Loop through provided number vector
  for (number in numbers) {
    
   # Create temporary df which stores result of 
   # applying function for multi metric probablity
   temp_df <- probability_metric_combos(df,
                              number,
                              n_sims) %>%
     mutate(no_par_left_out = (number - 1))
   
   # Append to list
   final_df <- c(final_df,
                 list(temp_df))
  }
  
  # Bind together dfs in the list
  final_df <- do.call(rbind,
                      final_df)
  
  return(final_df)
}

# How many do I want left out; 1 means all agree etc
range_of_metrics <- c(1,2,3)

# How many simulations per combo of noise, 
# gt and sample size
number_of_sims <- 250

# Calculate probability of agreement
agreement_prob <- leave_how_many_out(analysed_simulations,
                                     range_of_metrics,
                                     number_of_sims)

# Filter the data so that only 
# AICc, BIC and ANOVA remain
non_red_metrics_analysed <- analysed_simulations %>%
  dplyr::select(N_Obs,
                Ground_Truth_Model,
                SD_Noise,
                Sim_no,
                AICc_min,
                BIC_min,
                ANOVA_p_val_sel)

# Calculate probabilities for non redundant metrics
non_red_prob <- leave_how_many_out(
  non_red_metrics_analysed,
  c(1,2),
  250)

# Function for sorting when consensus 
# agrees on a particular number
metrics_agree_model <- function(df, 
                                n_met){
  
  # Filter out simulation number column 
  # as it is not needed
  df <- df %>%
    dplyr::select(-Sim_no)
  
  # Specify the columns to group by
  columns_to_group <- c("N_Obs", 
                        "Ground_Truth_Model", 
                        "SD_Noise")
  
  # Find all other columns
  other_cols <- setdiff(names(df), 
                        columns_to_group)
  
  # Function to calculate mode
  calculate_mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  
  # Group by N_Obs and SD_Noise and 
  # what the agreement between models is
  treat_df <- df %>%
    group_by(N_Obs, SD_Noise) %>%
    rowwise() %>%
    filter(n_distinct(c_across(other_cols)) <= n_met) %>%
    mutate(metric_consensus = calculate_mode(
      c_across(other_cols))) %>%
    group_by(metric_consensus)
  
  # df which contains rows where all 
  # metrics agree and are correct
  correct_df <- treat_df %>%
    filter(metric_consensus 
           == as.numeric(Ground_Truth_Model)) %>%
    ungroup() %>%
    dplyr::count(N_Obs,
                 metric_consensus,
                 SD_Noise) %>%
    dplyr::rename("n_corr" = "n")
  
  # count the columns for each mode
  agreement <- treat_df %>%
    ungroup() %>%
    dplyr::count(N_Obs,
                 metric_consensus,
                 SD_Noise)
  
  # df with only combos of models 
  # and noise to populate with probs
  only_vals <- df %>%
    group_by(N_Obs, SD_Noise) %>%
    rowwise() %>%
    filter(n_distinct(c_across(other_cols)) <= 3) %>%
    mutate(metric_consensus = calculate_mode(
      c_across(other_cols))) %>%
    ungroup() %>%
    dplyr::select(N_Obs,
                  SD_Noise,
                  metric_consensus) %>%
    unique()
  
  # combine dfs and calculate probs
  final_df <- list(only_vals,
                   agreement,
                   correct_df) %>%
    purrr::reduce(left_join,
                  by = c("N_Obs",
                         "metric_consensus",
                         "SD_Noise")) %>%
    mutate(prob_correct = ifelse(n != 0 & !is.na(n_corr),
                                 n_corr/n,
                                 0))
  
  return(final_df)}

# Apply so that all non redundant metrics agree
con_prob_1 <- metrics_agree_model(
  non_red_metrics_analysed,1) %>%
  ungroup() %>%
  mutate(no_par_left_out = 0)

# Apply so that all but one non redundant 
# metrics agree
con_prob_2 <- metrics_agree_model(
  non_red_metrics_analysed,2) %>%
  ungroup() %>%
  mutate(no_par_left_out = 1)

# Combine the dfs
consensus_prob <- rbind(con_prob_1,
                        con_prob_2)

# Filter for the n that you actually have for the study
prob_corr_final <- consensus_prob %>%
  filter(N_Obs == 10 
         & no_par_left_out == 1) %>%
  dplyr::select(-n,
                -n_corr,
                -no_par_left_out) %>%
  arrange(metric_consensus) %>%
  ungroup() %>%
  rename("Sample Size" = N_Obs,
         "Noise Level" = SD_Noise,
         "Selected Model" = metric_consensus,
         "Probability of Matching Ground Truth" = 
           prob_correct)

# Filter data to find deltaAIC values
delta_AICc_data <- sim_results[[3]] %>%
  group_by(N_Obs, 
           Ground_Truth_Model,
           SD_Noise,
           Sim_no) %>%
  mutate(delta_AICc_1_2 = ifelse(Model == "2",
                                 AICc - lag(AICc),
                                 NA)) %>%
  filter(Model == "2") %>%
  dplyr::select(delta_AICc_1_2)

delta_AICc_data_less <- delta_AICc_data %>%
  filter(abs(delta_AICc_1_2) < 2) %>%
  mutate(picks_model = ifelse(delta_AICc_1_2 > 0,
                          "1",
                          "2"),
         is_corr = ifelse(picks_model == 
                            Ground_Truth_Model,
                          T,
                          F))

delta_AICc_data_N <- delta_AICc_data_less %>%
  ungroup() %>%
  filter(N_Obs == 10) %>%
  group_by(N_Obs, SD_Noise, 
           Ground_Truth_Model) %>%
  summarise(
    prob_is_corr = mean(is_corr),
    n_total = n(),
    .groups = "drop"
  ) %>%
  mutate(prob_AICc_small = n_total/250)

#===================
# Reviewer 1 section
#===================
# Analyse the data generated using the ranges that R1
# requested
analysed_simulations_r1 <- analyse_simulations(
  sim_results_r1)

non_red_metrics_analysed_r1 <- analysed_simulations_r1 %>%
  dplyr::select(N_Obs,
                Ground_Truth_Model,
                SD_Noise,
                Sim_no,
                AICc_min,
                BIC_min,
                ANOVA_p_val_sel)

non_red_prob_r1 <- leave_how_many_out(
  non_red_metrics_analysed_r1,
  c(1,2),
  500)

con_prob_2_r1 <- metrics_agree_model(
  non_red_metrics_analysed_r1,
  2)

#============================
# Fractional CN sims analysis
#============================
# Apply the analysis functions
# Analyse fr CN sims w/og function
analysed_simulations_frCN <- analyse_simulations(
  sim_r_frCN)

# Run the truth finding function on
# frCN the simulated data
binary_truth_df_frCN <- calculate_individual_metric_truth(
  analysed_simulations_frCN)

# Probabilities on frCN data
prob_ind_metric_frCN <- calc_probabilities_individual_metrics(
  binary_truth_df_frCN)

# Calculate probability of agreement frCN
agreement_prob_frCN <- leave_how_many_out(
  analysed_simulations_frCN,
  range_of_metrics,
  100)

# Filter the data so that only 
# AICc, BIC and ANOVA remain
non_red_metrics_analysed_frCN <- analysed_simulations_frCN %>%
  dplyr::select(N_Obs,
                Ground_Truth_Model,
                SD_Noise,
                Sim_no,
                AICc_min,
                BIC_min,
                ANOVA_p_val_sel)

# Calculate probabilities for non redundant metrics
non_red_prob_frCN <- leave_how_many_out(
  non_red_metrics_analysed_frCN,
  c(1,2),
  250)

# Apply so that all non redundant metrics agree
con_prob_1_frCN <- metrics_agree_model(
  non_red_metrics_analysed_frCN,1) %>%
  ungroup() %>%
  mutate(no_par_left_out = 0)

# Apply so that all but one non redundant 
# metrics agree
con_prob_2_frCN <- metrics_agree_model(
  non_red_metrics_analysed_frCN,2) %>%
  ungroup() %>%
  mutate(no_par_left_out = 1)

# Combine the dfs
consensus_prob_frCN <- rbind(con_prob_1_frCN,
                        con_prob_2_frCN)

# Filter for the n that you actually have for the study
prob_corr_final_frCN <- consensus_prob_frCN %>%
  filter(N_Obs == 10 
         & no_par_left_out == 1) %>%
  dplyr::select(-n,
                -n_corr,
                -no_par_left_out) %>%
  arrange(metric_consensus) %>%
  ungroup() %>%
  rename("Sample Size" = N_Obs,
         "Noise Level" = SD_Noise,
         "Selected Model" = metric_consensus,
         "Probability of Matching Ground Truth" = 
           prob_correct)
