library(dplyr)
library(ggplot2)
library(ggrepel)
library(caret)
library(AICcmodavg)
library(patchwork)

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

#===========================
# Additional checks on model
#===========================
# Check residuals
residuals_lit_model <- ggplot(data = 
                                data.frame(
                                      
                                  comb_mod_lit$residuals),
                              aes(x = 
                                    comb_mod_lit.residuals)) +
  geom_histogram(binwidth = 1) +
  labs(y = "Count",
       x = "Literature Model Residual") +
  themething

# Mixed/random effects models for ligand
library(lme4)

# Filter for ligands that are represented more than 5 times
lit_dat_L_groups <- subsample_lit_df %>%
  filter(Ligand %in% c("PtBu3", 
                       "dppe", 
                       "RuPhos", 
                       "SIPr"))

#=============================================
# Models which cannot be accurately determined
# because of lack of data
#=============================================
# Mixed effects model for ligand
L_mix_eff_int_mod <- lmer(data = lit_dat_L_groups,
                        logk_rt ~ Soln_VIP_eV 
                        + V_bur 
                        + Nucleophilicity_eV 
                        + Electrophilicity_eV
                        + (1 | Ligand))

# Random effects model for ligand
L_uncor_mod <- lmer(data = lit_dat_L_groups,
                    logk_rt ~ Soln_VIP_eV 
                    + V_bur 
                    + Nucleophilicity_eV 
                    + Electrophilicity_eV 
                    + (Soln_VIP_eV | Ligand))

# Successful predictions for illustration
# Fixed effects for ligand aka varying intercept
L_fix_eff_mod <- lm(data = lit_dat_L_groups,
                    logk_rt ~ Soln_VIP_eV 
                    + V_bur 
                    + Nucleophilicity_eV 
                    + Electrophilicity_eV
                    + Ligand)

# Use VIP interaction to check for fixed intercept
# aka varying slope only
L_int_mod <- lm(data = lit_dat_L_groups,
                          logk_rt ~ Soln_VIP_eV 
                          + V_bur
                          + Nucleophilicity_eV 
                          + Electrophilicity_eV
                          + Ligand:Soln_VIP_eV)

# Add predictions from models to df
lit_dat_L_groups_pred <- lit_dat_L_groups %>%
  mutate(pred_logk_factor = predict(L_fix_eff_mod,
                                    lit_dat_L_groups),
         pred_res_factor = logk_rt - pred_logk_factor,
         pred_logk_int = predict(L_int_mod,
                                 lit_dat_L_groups),
         pred_res_int = logk_rt - pred_logk_int)

# Plot predicted v. observed
per_L_pred_plot <- lit_pred_plot(lit_dat_L_groups_pred %>%
                         mutate(pred_logk = 
                                  pred_logk_factor,
                                pred_res = pred_res_factor),
              -9,
              7,
              3) +
  scale_alpha(range = c(0.25,1),
              limits = c(-6.5,0))

#=======================================
# Reviewer 2 - Model comparison lit data
#=======================================
# Create the models
lit_E <- lm(data = subsample_lit_df,
            logk_rt ~ Soln_VIP_scaled
            + Nu_scaled + Electr_scaled)

lit_ES <- lm(data = subsample_lit_df,
            logk_rt ~ Soln_VIP_scaled 
            + Vbur_scaled 
            + Nu_scaled + Electr_scaled)

lit_ESCN <- lm(data = subsample_lit_df,
            logk_rt ~ Soln_VIP_scaled 
            + Vbur_scaled + Coord_no
            + Nu_scaled + Electr_scaled)

# Get the values for model selection
lit_mod_sel <- data.frame(
  model = c(1,2,3),
  BIC = c(BIC(lit_E),
          BIC(lit_ES),
          BIC(lit_ESCN)),
  AICc = c(AICc(lit_E),
           AICc(lit_ES),
           AICc(lit_ESCN)),
  pA = c(anova(lit_E,lit_ES)$`Pr(>F)`[[2]],
         anova(lit_ES,lit_ESCN)$`Pr(>F)`[[2]],
         anova(lit_E,lit_ESCN)$`Pr(>F)`[[2]])) %>%
  mutate(BIC_score = rank(-BIC),
         AICc_score = rank(-AICc),
         pA_score = c(1,3,2),
         model_sel_score = BIC_score 
         + AICc_score 
         + pA_score)

#============================
# Reviewer 2 VIP correlations
#============================
# Import Mulliken charge data and merge
lit_charge_data <- read.csv(
  "./literature_data/mulliken_charges_lit.csv") %>%
  mutate(ID = as.numeric(ID)) %>%
  left_join(subsample_lit_df,
            by = "ID")

# Plots for correlations
# HOMO and VIP
VIP_HOMO_lit_plot <- ggplot(data = NULL,
                            aes(x = Soln_VIP_eV,
                                y = Soln_E_HOMO)) +
  geom_smooth(data = lit_charge_data,
              method = "lm",
              color = "black") +
  geom_point(data = lit_charge_data %>%
               filter(Coord_no == 4),
             fill = "#B71B1B",
             size = 3,
             pch = 21) +
  geom_point(data = lit_charge_data %>%
               filter(Coord_no == 3),
             fill = "#636C9D",
             size = 3,
             pch = 21) +
  labs(x = "VIP in Solution (eV)",
       y = "HOMO Energy in Solution (Hartrees)") +
  themething

# Mulliken charges
VIP_mulliken_lit <- ggplot(data = NULL,
                           aes(x = Soln_VIP_eV,
                               y = SPVIP_Charge_Block1)) +
  geom_smooth(data = lit_charge_data,
              method = "lm",
              color = "black") +
  geom_point(data = lit_charge_data %>%
               filter(Coord_no == 4),
             fill = "#B71B1B",
             size = 3,
             pch = 21) +
  geom_point(data = lit_charge_data %>%
               filter(Coord_no == 3),
             fill = "#636C9D",
             size = 3,
             pch = 21) +
  labs(x = "VIP in Solution (eV)",
       y = "Mulliken Charge") +
  themething

# Combine the plots
VIP_corr_plots <- VIP_mulliken_lit +
  VIP_HOMO_lit_plot +
  plot_layout(axes = "collect")

print(VIP_corr_plots)

# VIP/HOMO relationship
VIP_HOMO_lm <- lm(data = lit_charge_data,
                  Soln_E_HOMO ~ Soln_VIP_eV)

# Extract the residuals
lit_res_df <- lit_charge_data %>%
  na.omit() %>%
  mutate(HOMO_VIP_residuals = VIP_HOMO_lm$residuals)

# Plot the residuals vs VIP
HOMO_VIP_res_plot <- ggplot(data = NULL,
                            aes(x = Soln_VIP_eV,
                                y = 
                                  HOMO_VIP_residuals)) +
  geom_abline(intercept = 0,
              slope = 0) +
  geom_point(data = lit_res_df %>%
               filter(Coord_no == 4),
             fill = "#B71B1B",
             size = 3,
             pch = 21) +
  geom_point(data = lit_res_df %>%
               filter(Coord_no == 3),
             fill = "#636C9D",
             size = 3,
             pch = 21) +
  labs(x = "VIP in Solution (eV)",
       y = "Residual with HOMO") +
  themething

# Histogram of residuals
HOMO_VIP_res_hist <- ggplot(data = lit_res_df) +
  geom_histogram(aes(x = HOMO_VIP_residuals),
                 binwidth = 0.001) +
  labs(y = "Count",
       x = "HOMO/VIP Residual") +
  themething

# Combine residual plots
res_plots <- HOMO_VIP_res_hist +
  HOMO_VIP_res_plot

# Look for outliers
outlier_lit_res_df <- lit_res_df %>%
  merge(subsample_lit_df %>%
              mutate(
    model_cd  = cooks.distance(comb_mod_lit_scaled),
    model_mean_cd = mean(model_cd))) %>%
  mutate(HOMO_VIP_cd = cooks.distance(VIP_HOMO_lm),
         HOMO_VIP_mean_cd = mean(HOMO_VIP_cd))

# Plot the two cook's distances
cd_plot <- ggplot(data = NULL,
                  aes(y = HOMO_VIP_cd,
                      x = model_cd)) +
  geom_point(data = outlier_lit_res_df %>%
               filter(Coord_no == 4),
             fill = "#B71B1B",
             size = 3,
             pch = 21) +
  geom_point(data = outlier_lit_res_df %>%
               filter(Coord_no == 3),
             fill = "#636C9D",
             size = 3,
             pch = 21) +
  labs(x = "Model Cook's Distance",
       y = "VIP/HOMO Cook's Distance") +
  themething

# Model residuals
res_lit_sample_df <- subsample_lit_df %>%
  mutate(model_residuals = 
           comb_mod_lit_scaled$residuals)

# Plot model residuals vs log10k
res_logk_plot <- ggplot(data = NULL,
                  aes(y = model_residuals,
                      x = logk_rt)) +
  geom_abline(intercept = 0,
              slope = 0) +
  geom_point(data = res_lit_sample_df %>%
               filter(Coord_no == 4),
             fill = "#B71B1B",
             size = 3,
             pch = 21) +
  geom_point(data = res_lit_sample_df %>%
               filter(Coord_no == 3),
             fill = "#636C9D",
             size = 3,
             pch = 21) +
  labs(x = expression(paste("Measured ",
                            log[10](k),
                            " at 25 \u00B0C")),
       y = "Literature Model Residual") +
  themething

model_residuals_plots <- residuals_lit_model +
  res_logk_plot
