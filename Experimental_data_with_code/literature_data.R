library(dplyr)
library(ggplot2)
library(ggrepel)
library(caret)

#============
# Import data
#============
subsample_lit_df <- read.csv(
  file.path(
    "./literature_data/subsample_lit.csv")) %>%
  mutate(T_kelvin = Temperature_C + 273.15,
         Coord_no = as.factor(Coord_no),
         k_rt = if_else(T_kelvin != 298.15, 
                        2.1e10*
                          T_kelvin*
                          exp((T_kelvin/298.15)*
                                log(
                          Rate_constant_inverse_seconds/
                                      (2.1e10*T_kelvin))), 
                        Rate_constant_inverse_seconds),
         logk_rt = log10(k_rt),
         Soln_VIP_scaled = scale(Soln_VIP_eV,
                                 center = T,
                                 scale = T),
         Gas_VIP_scaled = scale(Gas_VIP_eV,
                                 center = T,
                                 scale = T),
         Vbur_scaled = scale(V_bur,
                             center = T,
                             scale = T),
         sol_E_HOMO_scaled = scale(Soln_E_HOMO,
                                center = T,
                                scale = T),
         Nu_scaled = scale(Nucleophilicity_eV,
                           center = T,
                           scale = T),
         Electr_scaled = scale(Electrophilicity_eV,
                               center = T,
                               scale = T))

#=============
# Apply models
#=============
# Raw model
comb_mod_lit <- lm(data = subsample_lit_df,
                       logk_rt ~ Soln_VIP_eV + V_bur 
                       + Nucleophilicity_eV 
                   + Electrophilicity_eV)

# Scaled model
comb_mod_lit_scaled <- lm(data = subsample_lit_df,
                   logk_rt ~ Soln_VIP_scaled 
                   + Vbur_scaled 
                   + Nu_scaled + Electr_scaled)

# Create df w/predictions from model
lit_pred_df <- subsample_lit_df %>%
  mutate(pred_logk = predict(comb_mod_lit_scaled,
                             subsample_lit_df),
         pred_res = logk_rt - pred_logk)

# Plot of predicted vs. observed
pred_plot_meta_regression <- ggplot(data = lit_pred_df,
                                    aes(x = logk_rt,
                                        y = pred_logk,
                                        shape = Bond_formed,
                                        alpha = 
                                          -abs(pred_res),
                                        fill = Coord_no)) +
  geom_abline(slope = 1,
              intercept = 0,
              color = "#4E525A",
              linewidth = 1) +
  geom_abline(slope = 1,
              intercept = -1,
              color = "#4E525A",
              linewidth = 0.5,
              alpha = 0.5,
              linetype = "dashed") +
  geom_abline(slope = 1,
              intercept = 1,
              color = "#4E525A",
              linewidth = 0.5,
              alpha = 0.5,
              linetype = "dashed") +
  geom_point(size = 3) +
  scale_fill_manual(values = c("#636C9D",
                               "#B71B1B")) +
  scale_shape_manual(values = c(21,22,23)) +
  labs(x = expression(paste("Measured ",
                            log[10](k),
                            " at 25 \u00B0C")),
       y = expression(paste("Predicted ", 
                            log[10](k),
                            " at 25 \u00B0C")),
       fill = "Coord. No.",
       shape = "Bond Formed",
       alpha = "Residisual") +
  guides(fill = guide_legend(override.aes = 
                               list(color = c("#636C9D",
                                             "#B71B1B")))) +
  themething

#=================
# Model validation
#=================
# Test train split w/caret package
tt_lit <- train(logk_rt ~ Soln_VIP_scaled 
                + Vbur_scaled 
                + Nu_scaled + Electr_scaled,
      method = "lm",
      data = lit_pred_df,
      trControl = trainControl(
        method = "LGOCV",
        p = 0.6,
        number = 50,
        savePredictions = "final"))

# Extract the predictions from the cross-validation
pred_train <- tt_lit$pred

# Filter the test predictions
test_pred <- 
  pred_train[pred_train$Resample 
              != pred_train$Resample[1], ]

# Calculate the R-squared for each resample
r2_values <- by(test_pred, 
                test_pred$Resample, 
                function(df) {
                  cor(df$obs, df$pred)^2
                })

# Display the mean test R-squared value across 
# all resamples
mean_r2_test <- mean(unlist(r2_values))

# LOOCV function
loocv_lit <- train(logk_rt ~ Soln_VIP_scaled 
                   + Vbur_scaled 
                   + Nu_scaled + Electr_scaled,
                   method = "lm",
                      data = lit_pred_df,
                      trControl = trainControl(
                        method = "LOOCV"))

#===========
# Outlier ID
#===========
# Create a plot
lit_pred_add_cd <- lit_pred_df %>%
  mutate(Cooks_d  = cooks.distance(comb_mod_lit_scaled),
         mean_c_d = mean(Cooks_d))

# Plot Cook's D
cd_plot_lit <- ggplot(data = lit_pred_add_cd %>%
                        filter(Cooks_d > 3*mean_c_d)) +
  geom_bar(aes(x = as.factor(ID),
               y = Cooks_d),
           stat = "identity",
           fill = "red") +
  geom_hline(yintercept = 3*lit_pred_add_cd$mean_c_d,
             linetype = "dashed",
             alpha = 0.3) +
  geom_hline(yintercept = 4/(nrow(lit_pred_add_cd)
                             - 4 - 1),
             linetype = "dashed",
             alpha = 0.3,
             color = "red") +
  annotate("text",
           x = "10",
           y = lit_pred_add_cd$mean_c_d*2.5,
           label = "Conservative \n Outlier \n Threshdold",
           size = 4.5,
           alpha = 0.3) +
  annotate("text",
           x = "10",
           y = (4/(nrow(lit_pred_add_cd)
                  - 4 - 1))*1.3,
           label = "Outlier \n Threshdold",
           size = 4.5,
           alpha = 0.3,
           color = "red") +
  labs(y = "Cook's Distance",
       x = "ID") +
  themething

#=========================
# Re-evaluate w/o outliers
#=========================
# Remove outliers
lit_df_filt_outliers <- lit_pred_add_cd %>%
  filter(Cooks_d < (4/(nrow(lit_pred_add_cd) - 4 - 1)))

# Test train split w/caret package no outliers
tt_lit_no_out <- train(logk_rt ~ Soln_VIP_scaled + 
                         Vbur_scaled + Nu_scaled + 
                         Electr_scaled,
                method = "lm",
                data = lit_df_filt_outliers,
                trControl = trainControl(
                  method = "LGOCV",
                  p = 0.6,
                  number = 50,
                  savePredictions = "final"))

# Apply LOOCV to data w/o outliers
loocv_lit_no_out <- train(logk_rt ~ Soln_VIP_scaled + 
                            Vbur_scaled + Nu_scaled + 
                            Electr_scaled,
                          method = "lm",
                            data = lit_df_filt_outliers,
                            trControl = trainControl(
                              method = "LOOCV"))

# Extract the predictions from the cross-validation
predictions <- tt_lit_no_out$pred

# Filter the test predictions
test_predictions <- 
  predictions[predictions$Resample != predictions$Resample[1], ]

# Calculate the R-squared for each resample
r2_values <- by(test_predictions, test_predictions$Resample, 
                function(df) {
  cor(df$obs, df$pred)^2
})

# Display the mean test R-squared value 
# across all resamples
mean_r2 <- mean(unlist(r2_values))

#========================
# Publication-level plots
#========================
# Fucntion for making a predicted plot
lit_pred_plot <- function(inp_df,
                          low_lim,
                          high_lim,
                          pt_size){
  
  out_plot <- ggplot(data = inp_df,
                     aes(x = logk_rt,
                         y = pred_logk,
                         shape = Bond_formed,
                         alpha = 
                           -abs(pred_res),
                         fill = Coord_no)) +
    geom_abline(slope = 1,
                intercept = 0,
                color = "#4E525A",
                linewidth = 1) +
    geom_abline(slope = 1,
                intercept = -1,
                color = "#4E525A",
                linewidth = 0.5,
                alpha = 0.5,
                linetype = "dashed") +
    geom_abline(slope = 1,
                intercept = 1,
                color = "#4E525A",
                linewidth = 0.5,
                alpha = 0.5,
                linetype = "dashed") +
    geom_point(size = pt_size) +
    xlim(low_lim,high_lim) +
    ylim(low_lim,high_lim) +
    labs(x = expression(paste("Measured ",
                              log[10](k),
                              " at 25 \u00B0C")),
         y = expression(paste("Predicted ", 
                              log[10](k),
                              " at 25 \u00B0C")),
         fill = "Coord. No.",
         shape = "Bond Formed",
         alpha = "Pred. Resid.") +
    scale_fill_manual(values = c("#636C9D",
                                  "#B71B1B")) +
    scale_shape_manual(values = c(21,22,23)) +
    themething
  
  return(out_plot)
}

# Prediction plot
lit_pred_plot_pub <- lit_pred_plot(lit_df_filt_outliers,
                                   -9,
                                   7,
                                   2) +
scale_alpha(range = c(0.25,1),
            limits = c(-6.5,0))

# By ligand
lig_grid_plot <- lit_pred_plot(lit_pred_df %>%
  filter(Ligand %in% c("PtBu3",
                       "SIPr",
                       "RuPhos",
                       "dppe")),
                               -11,8,3) +
  facet_wrap(~ Ligand) +
  theme(legend.position = "none")

# By ligand for manuscript figure
lig_grid_plot_ms <- lit_pred_plot(lit_pred_df %>%
    filter(Ligand %in% c("PtBu3",
                         "dppe")),
                               -8,7,3) +
  facet_wrap(~ Ligand) +
  theme(legend.position = "none") +
  scale_alpha(range = c(0.25,1),
              limits = c(-5,0))

# Apply the model to each group of ligand 
# that has enough data
# Make function to do so
apply_lit_mod_per_L <- function(inp_df,
                                L){
  
  df_filt <- inp_df %>%
    filter(Ligand == as.character({{L}}))
  
  mod_filt <- lm(data = df_filt,
                               logk_rt ~ Soln_VIP_scaled 
                               + Vbur_scaled 
                               + Nu_scaled + Electr_scaled)
  
  df_filt_pred <- df_filt %>%
    mutate(pred_logk = predict(mod_filt,
                                 df_filt))
  
  out_list <- list(df_filt_pred,
                   summary(mod_filt))
  
  return(out_list)
}

# PtBu3
lit_filt_PtBu3 <- apply_lit_mod_per_L(subsample_lit_df,
                                      "PtBu3")

# SIPr
lit_filt_SIPr <- apply_lit_mod_per_L(subsample_lit_df,
                                      "SIPr")

# dppe
lit_filt_dppe <- apply_lit_mod_per_L(subsample_lit_df,
                                     "dppe")

# RuPhos
lit_filt_RuPhos <- apply_lit_mod_per_L(subsample_lit_df,
                                     "RuPhos")

# Combine data to plot
lit_pred_df_per_L <- lit_filt_PtBu3[[1]] %>%
  rbind(lit_filt_SIPr[[1]],
        lit_filt_dppe[[1]],
        lit_filt_RuPhos[[1]]) %>%
  mutate(pred_res = pred_logk - logk_rt)

# Create plot
per_L_adj_mod_plot <- lit_pred_plot(lit_pred_df_per_L,
                                    -9,
                                    6,
                                    3) +
  facet_wrap(~ Ligand) +
  theme(legend.position = "none")

# Only two Ls for ms figure
per_L_adj_mod_plot_ms <- lit_pred_plot(
  lit_pred_df_per_L %>%
    filter(Ligand %in%
             c("dppe",
               "PtBu3")),
                                    -8,
                                    7,
                                    3) +
  facet_wrap(~ Ligand) +
  theme(legend.position = "none") +
  scale_alpha(range = c(0.25,1),
              limits = c(-5,0))