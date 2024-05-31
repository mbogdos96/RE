library(MASS)
library(caret)
library(dplyr)
library(AICcmodavg)

#================
# Data generation
#================

generate_data <- function(n, noise_sd, ground_truth_model, seed) {
  set.seed(as.numeric(seed))  # for reproducibility
  
  # Sample oxidation potentials between 0 and 2
  E_ox <- runif(n, 0.1, 2)
  
  # Sample buried volumes between 20 and 80
  V_bur <- runif(n, 20, 80)
  
  # Ensure at least two instances of each level of L_type, 
  #then fill the rest as close to 1:1 as possible
  if(n > 2) {
    L_type <- factor(c("mono", "mono", 
                       "bi", "bi", 
                       sample(c("mono", "bi"), 
                              n-4, replace = T)))
  } else {
    # For cases where n <= 2, just alternate between mono and bi
    L_type <- factor(rep(c("mono", "bi"), 
                         length.out = n))
  }
  
  # Create a df with the data, scale and center
  data <- data.frame(E_ox, 
                     V_bur, 
                     L_type) %>%
    mutate(E_ox_scaled = scale(E_ox,
                               center = T,
                               scale = T),
           V_bur_scaled = scale(V_bur,
                                center = T,
                                scale = T))
  
  # Generate noise of a given amount
  epsilon <- rnorm(n, 
                   0, 
                   noise_sd)
  
  # Generate random coefficient for intercept
  intercept <- runif(1, -10, -0.1)
  coef_E_ox <- runif(1, 1, 10)
  coef_V_bur <- runif(1, -10, -1)
  coef_L_type <- runif(1, 0.1, 10)
  
  # Generate the logk data based on different ground truths
  if (ground_truth_model == 1) {
    logk <- intercept + 
      coef_E_ox * data$E_ox_scaled + 
      epsilon
    
    ground_truth <- c(intercept = intercept, 
                      coef_E_ox = coef_E_ox)
  } else if (ground_truth_model == 2) {
    logk <- intercept + 
      coef_E_ox * data$E_ox_scaled + 
      coef_V_bur * data$V_bur_scaled + 
      epsilon
    ground_truth <- c(intercept = intercept, 
                      coef_E_ox = coef_E_ox, 
                      coef_V_bur = coef_V_bur)
  } else if (ground_truth_model == 3) {
    logk <- intercept + 
      coef_E_ox * data$E_ox_scaled + 
      coef_V_bur * data$V_bur_scaled + 
      coef_L_type * ifelse(data$L_type == "mono", 
                           1, 
                           0) + 
      epsilon
    ground_truth <- c(intercept = intercept, 
                      coef_E_ox = coef_E_ox, 
                      coef_V_bur = coef_V_bur, 
                      coef_L_type = coef_L_type)
  }
  
  # Save logk into df
  data_df <- data.frame(logk) %>%
    cbind(data)
  
  # Outputs
  list_out <- list(data_df,
                   ground_truth)
  
  # Output dfs
  return(list_out)
}

#=====================================
# Model fitting and metric calculation
#=====================================

fit_models <- function(input) {
  data <- input[[1]]
  ground_truth <- input[[2]]
  
  # Ensure L_type is treated as factor with both levels present
  data$L_type <- factor(data$L_type, 
                        levels = c("mono", "bi"))
  
  # Explicitly define contrasts for L_type
  contrasts(data$L_type) <- contr.treatment(2)
  
  # Define the models to be compared
  model1 <- lm(logk ~ E_ox_scaled, 
               data = data)
  model2 <- lm(logk ~ E_ox_scaled + V_bur_scaled,
               data = data)
  model3 <- lm(logk ~ E_ox_scaled + V_bur_scaled + L_type, 
               data = data)
  
  # Create data frame for metrics
  metrics <- data.frame(
    Model = c("1", 
              "2", 
              "3"),
    AICc = c(AICc(model1), 
             AICc(model2), 
             AICc(model3)),
    BIC = c(BIC(model1), 
            BIC(model2), 
            BIC(model3)),
    adj_R2 = c(summary(model1)$adj.r.squared, 
               summary(model2)$adj.r.squared, 
               summary(model3)$adj.r.squared),
    Mean_Prediction_Interval_Width = c(mean((
      predict(model1, 
              interval = "prediction")[,3]
      - predict(model1,
                interval = "prediction")[,2])),
      mean((predict(model2, 
                    interval = "prediction")[,3]
            - predict(model2,
                      interval = "prediction")[,2])),
      mean((predict(model3,
                    interval = "prediction")[,3]
            - predict(model3,
                      interval = "prediction")[,2]))),
    ANOVA_log10p = c(log10(anova(model1, model2)$Pr[2]), 
                     log10(anova(model2, model3)$Pr[2]), 
                     log10(anova(model1, model3)$Pr[2])),
    Estimate_E_ox = c(summary(model1)$coef[2, "Estimate"], 
                      summary(model2)$coef[2, "Estimate"], 
                      summary(model3)$coef[2, "Estimate"]),
    SE_Estimate_E_ox = c(summary(model1)$coef[2, "Std. Error"], 
                          summary(model2)$coef[2, "Std. Error"], 
                          summary(model3)$coef[2, "Std. Error"]),
    Estimate_Intercept = c(summary(model1)$coef[1, "Estimate"], 
                           summary(model2)$coef[1, "Estimate"], 
                           summary(model3)$coef[1, "Estimate"]),
    SE_Estimate_Intercept = c(summary(model1)$coef[1, "Std. Error"], 
                          summary(model2)$coef[1, "Std. Error"], 
                          summary(model3)$coef[1, "Std. Error"]))
  
  # Add ground truth coefficients
  for (coef_name in names(ground_truth)) {
    metrics[[paste0("GT_", coef_name)]] <- rep(ground_truth[[coef_name]], 3)
  }

  metrics <- metrics %>%
    mutate(Delta_coef_GT_E_ox = (Estimate_E_ox - GT_coef_E_ox),
           Delta_coef_GT_Intercept = (Estimate_Intercept - GT_intercept))
  
  return(metrics = metrics)
}

#================
# Run simulations
#================

run_simulations <- function(n_obs_vec, noise_sd_vec, n_simulations) {
  
  # Initialize an empty dataframe to store simulated data
  sim_data_list <- list()  
  
  # Initialize df to store the metrics for the datasets
  metrics_list <- list()
  
  # Initialize an empty list to store seeds
  seeds_list <- list()
  
  for (n_obs in n_obs_vec) {
    # Iterate over noise_sd values
    for (noise_sd in noise_sd_vec) {  
      for (ground_truth_model in 1:3) {
        for (i in 1:n_simulations) {
          # Generate a single seed for this combination
          seed <- sample(10000000, 1)  
          
          # Store seed in the list
          seeds_list <- c(seeds_list, 
                          list(seed))
          
          # Initialize dataframe to store metrics for this combination
          metrics_df <- data.frame()
          
          # Set seed for reproducibility
          set.seed(seed)  
          
          # Generate the data using the function above
          data_sim <- generate_data(n_obs, 
                                noise_sd, 
                                ground_truth_model, 
                                seed)
          
          # Add other information to the dataframe
          data_sim[[1]]$SD_Noise <- noise_sd
          data_sim[[1]]$Ground_Truth_Model <- ground_truth_model
          data_sim[[1]]$N_Obs <- n_obs
          data_sim[[1]]$Sim_no <- i
          
          # Append simulated data to df
          sim_data_list <- c(sim_data_list,
                             list(data_sim[[1]]))
          
          # Calculate the metrics for the simulated data
          metrics <- fit_models(data_sim)
          
          # Add other information to the dataframe
          metrics$SD_Noise <- noise_sd
          metrics$Ground_Truth_Model <- ground_truth_model
          metrics$N_Obs <- n_obs
          metrics$Sim_no <- i
        
          # Append metrics to results dataframe
          metrics_list <- c(metrics_list,
                            list(metrics))
        }
    }
  }
 }
  
  sim_data <- bind_rows(sim_data_list)
  
  results <- bind_rows(metrics_list)
  
  return(list(simulated_data = sim_data, 
              seeds_list = seeds_list, 
              metrics = results))
}

# Parameters to simulate across
n_obs_vec <- seq(5, 30, by = 1)
noise_sd <- c(0.25,1,2)
n_simulations <- 250

# Run simulations
sim_results <- run_simulations(n_obs_vec, 
                               noise_sd, 
                               n_simulations)