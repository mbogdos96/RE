library(dplyr)
library(ggplot2)
library(nls2)
library(ggrepel)
library(propagate)
library(jtools)
library(patchwork)
library(sjPlot)
library(modelsummary)

source("SW_plots.R")
source("common_functions.R")

# Define the molecular weight of the 
# internal standard used in all cases
internal_standard_molecular_weight <- 246.22

#==============
#==============
# Data handling
#==============
#==============

#====================================
# Charton-Taft parameter calculations
#====================================
# Import data
taft_ref_data <- read.csv(file.path(
  "./half_life_data/taft_19_f.csv"))

# 19F delta delta and sigma i regression
lm_19F_sigma <- lm(data = taft_ref_data,
                   sigma_i ~ delta_F)

# Enter delta delta 19F data predicted for arenes
F4_aryl <- data.frame(delta_F = c(0.74,
                                  -2.63)) 

# Predict values based on regression
F4_aryl_pred <- F4_aryl %>%
  mutate(sigma_i = predict(lm_19F_sigma,
                           F4_aryl)) 

# Sigma i for CF3 arene using COOH regression method
# Import data collected from literature
taft_pKa_data <- read.csv(file.path(
  "./half_life_data/taft_pKa.csv"))

# Regression of sigma i on pKa in water from lit
lm_pKa_sigma <- lm(data = taft_pKa_data,
                   sigma_i ~ pKa_H2O)

# Enter the predicted pKa for the relevant arene
bisCF3_aryl <- data.frame(pKa_H2O = c(3.34)) 

# Predict the value sigma i CF3 arene based on
# regression
bisCF3_aryl_pred <- bisCF3_aryl %>%
  mutate(sigma_i = predict(lm_pKa_sigma,
                           bisCF3_aryl))

# Sigma i for ClMs using interpolation
# Import literature data
sigma_sulf <- read.csv(file.path(
  "./half_life_data/taft_sulf_alk.csv"))

# Regression between two sigma i values
lm_sulf_sigma_i <- lm(data = sigma_sulf,
                      sigma_i_sulf ~ sigma_i_alk)

# Enter data reported for ClMs alkyl sigma i
clMs_sigma_i <- data.frame(sigma_i_alk = 0.13) 

# Predict value based on regression
clMs_sigma_i_pred <- clMs_sigma_i %>%
  mutate(sigma_i_sulf = predict(lm_sulf_sigma_i,
                                clMs_sigma_i))

#=======================================
# Calculate half-lives for palladacycles
#=======================================
# Import the data and then convert values to SI units
half_life_data <- read.csv(
  file.path(
    "./half_life_data/half_lives_backbones_PtBu3.csv")) %>%
                  mutate(
      Temperature_K = Temperature_C + 273.15,
      time_s = Time_min*60,
      Internal_standard_mol = 
      (Internal_standard_mass_mg*0.001)/
      internal_standard_molecular_weight,
      Pd_II_mol= Pd_II_integral*Internal_standard_mol) %>%
      group_by(Backbone) %>%
      mutate(Pd_II_mol_norm = 
               Pd_II_mol/max(Pd_II_mol)) %>%
      ungroup()

# Function which fits first order 
# rate equation to estimate a rate constant
calculate_rate_constant <- function(df,
                                    x,
                                    y) {
  # Get the names of backbones in a list
  unique_backbones <- unique(df$Backbone)
  
  # Empty df to store results
  result <- data.frame()
  
  # Loop through backbone list
  for (backbone in unique_backbones) {
    
    # Keep only data for each backcbone in 
    # input df
    subset_data <- filter(df, 
                          Backbone == backbone)
    
    # Fit equation to data for that backbone
    model_fit <- try(
      nls(Pd_II_mol ~ A * exp(-B * time_s), 
          data = subset_data, 
          start = list(A = x, 
                       B = y)
      )
    )
    
    # Loop for handling errors thrown if
    # provided start values for coefs
    # lead to issues
    if (!inherits(model_fit, 
                  "try-error")) {
      # If if succeeds extract the coef values
      coefficients <- coef(model_fit)
      
      # Store these in result df
      result <- bind_rows(result, 
                          data.frame(
                            Backbone = backbone,
                            A0 = coefficients[["A"]], 
                            rate_constant = 
                              coefficients[["B"]]))
    }
  }
  
  # Append results
  result <- result
  
  # Return results
  return(result)
}

# Apply the function to the backbone data, 
# manipulate df and calculate t half at rt
backbone_half_lives_PtBu3 <- calculate_rate_constant(
  half_life_data,
  0.1,
  1e-7) %>%
  left_join(dplyr::select(half_life_data,
                          Backbone,
                          N_H_stretch_IR,
                          sigma_para,
                          Temperature_K,
                          sigma_i_n,
                          sigma_i_ar), 
            by = "Backbone") %>%
  distinct() %>%
  mutate(rate_constant = if_else(Temperature_K != 298.15, 
                                 2.1e10*
                                   Temperature_K*
                      exp((Temperature_K/298.15)*
                               log(rate_constant/
                            (2.1e10*Temperature_K))), 
                                 rate_constant),
         half_life_s = log(2)/rate_constant,
         log_half_life = log10(half_life_s),
         sum_s_i_a = sigma_i_n + sigma_i_ar,
         sum_s_i_a_scaled = sum_s_i_a,
         sigma_para_scaled = as.numeric(sigma_para),
         sigma_i_n_scaled = as.numeric(sigma_i_n),
         sigma_i_ar_scaled = as.numeric(sigma_i_ar),
         NH_IR_scaled = as.numeric(N_H_stretch_IR),
         across(c(sigma_para_scaled, 
                  sigma_i_n_scaled,
                  sigma_i_ar_scaled,
                  NH_IR_scaled,
                  sum_s_i_a_scaled),
                ~ scale(.,
                        scale = T,
                        center = T)))

# Remove NH from data before building model, 
# as it is outwith measurement range for rt NMR
t_half_model <- backbone_half_lives_PtBu3 %>%
  mutate(logk = log10(rate_constant),
         Backbone = na_if(Backbone,"CF3-NH-Ph"),
         Special_label_exp = sapply(
           Backbone,
           replace_backbone_labels)) %>%
  na.omit()

# logt_half NH IR and sigma p model
lm_backbones_IR <- lm(data = t_half_model,
  log_half_life ~  sigma_para_scaled + NH_IR_scaled)

# logt_half sigma model 3 param
lm_backbones_sigma <- lm(data = t_half_model,
 log_half_life ~  sigma_para_scaled +
   sigma_i_ar_scaled +
   sigma_i_n_scaled)

# logt_half sigma model 2 param
lm_sum_back <- lm(data = t_half_model,
 log_half_life ~  sigma_para_scaled +
   sum_s_i_a_scaled)

# Predict w/sum model
t_half_pred_df_2 <- t_half_model %>%
  mutate(log_t_half_pred_sum = predict(lm_sum_back,
                                       t_half_model %>%
                                         filter(
                                Backbone != "CF3-NH-Ph"
                                         )))

#========================
# Alternative descriptors
#========================
# Import data
nu_e_param_fragments <- read.csv(file.path(
  "./half_life_data/nu_e_param.csv"
))

# Function for extracting the right values
extr_e_nu_params <- function(df_params,
                             df_append){
  # Add new columns for E_HOMO_amido and aryl_electr
  df_append$E_HOMO_amido <- NA
  df_append$aryl_electr <- NA
  
  # Loop through each Backbone in sum_test_df
  for (i in 1:nrow(df_append)) {
    # Split the backbone into fragments
    fragments <- unlist(strsplit(df_append$Backbone[i], 
                                 "-"))
    
    # Amido fragment join by '-'
    amido_fragment <- paste(fragments[1], 
                            fragments[2], 
                            sep = "-")
    
    # Aryl fragment (last element)
    aryl_fragment <- fragments[3]
    
    # Find E.HOMO...eV for amido_fragment
    amido_match <- df_params$Fragment == amido_fragment
    if (any(amido_match)) {
      df_append$E_HOMO_amido[i] <- 
        df_params$E.HOMO...eV[amido_match]
    }
    
    # Find electrophilicity...eV for aryl_fragment
    aryl_match <- df_params$Fragment == aryl_fragment
    if (any(aryl_match)) {
      df_append$aryl_electr[i] <- 
        df_params$electrophilicity...eV[aryl_match]
    }
  }
  
  return(df_append)
}

# Apply to half-life df
hl_e_nu_df <- extr_e_nu_params(nu_e_param_fragments,
                               t_half_pred_df_2) %>%
  filter(Backbone != "CF3-NH-Ph") %>%
  mutate(Nu_scaled = scale(as.numeric(E_HOMO_amido),
                           scale = T,
                           center = T),
         Electr_scaled = scale(as.numeric(aryl_electr),
                               scale = T,
                               center = T),
         Special_label_exp = sapply(
           Backbone,
           replace_backbone_labels))

# Linear model w/electrophilicity and E HOMO
e_nu_lm_backbones <- lm(data = hl_e_nu_df,
                        log_half_life ~ Nu_scaled + 
                          Electr_scaled)

#=========
#=========
# Plotting
#=========
#=========

#=========
# Kinetics
#=========
# Plot showing concentrations over time 
# for the various backbones
backbone_t_plot <- ggplot(data = half_life_data,
                          aes(x = time_s,
                              y = Pd_II_mol_norm)) +
  geom_point(aes(color = Backbone),
             size = 3) +
  scale_color_brewer(palette = "Paired") +
  labs(x = "Time (s)",
       y = "Relative Pd(II) Concentration") +
  themething

#==============
# Distributions
#==============
# Histogram of half lives
t_half_hist <- ggplot(data = 
                        backbone_half_lives_PtBu3,
                      aes(x = half_life_s)) +
  geom_histogram(binwidth = 1e6) +
  geom_density(fill = "orange",
               alpha = 0.2) +
  labs(x = "Half Life at 25 \u00B0C (s)",
       y = "Count") +
  themething

# Histogram log half lives
log_t_half_hist <- ggplot(data = 
                            backbone_half_lives_PtBu3,
                          aes(x = log_half_life)) +
  geom_histogram(binwidth = 1) +
  geom_density(fill = "orange",
               alpha = 0.2) +
  labs(x = expression(paste(log[10]("Half-life"), 
                            " at 25 \u00B0C")),
       y = "Count") +
  themething

# Combine the two histograms into one plot
hist_comb <- t_half_hist +
  log_t_half_hist +
  plot_layout(axis_titles = "collect")

#=========================
# Charton-Taft regressions
#=========================
# Function for plotting the correlations 
# from which I calculate the sigma_i
# values for the various substituents - needed for SI
# The arguments it takes are 
# df_data - data frame your data is in
# x_var - whatever you want plotted on the x axis; 
# y_var - cf. x_var
# pt_size - size of points, pt_col - color of points, 
# x_lab - string 
# that is x axis label; y_lab - cf. x_lab
simple_plot_mb_style <- function(df_data,
                                 x_var,
                                 y_var,
                                 pt_size,
                                 pt_col,
                                 x_lab,
                                 y_lab){
  plot <- ggplot(data = df_data,
                 aes(x = {{x_var}},
                     y = {{y_var}})) +
    geom_smooth(method = lm,
                se = T,
                inherit.aes = T) +
    geom_point(size = pt_size,
               color = pt_col) +
    labs(x = x_lab,
         y = y_lab) +
    themething
}

# Plot for correlation between 19F NMR signals 
# from Taft paper
taft_19F_plot <- simple_plot_mb_style(taft_ref_data,
                                      delta_F,
                                      sigma_i,
                                      3,
                                      "black",
                                      expression(
                                        paste(
                                          {}^19,
                                          "F{",
                                          {}^1,
              "H} NMR signal relative to PhF")),
                                      expression(
                                        sigma[I]))

# Regression table for Taft 19F
Taft_ref_table <- modelsummary(lm_19F_sigma,
  fmt = 3,
  estimate = "{estimate} [{conf.low}, {conf.high}]",
  statistic = NULL,
  gof_omit = "AIC|BIC|Log.Lik.",
  title = "Summary Statistics for Taft 19F Model",
  output = "latex")

# Plot for correlation between pKa and sigma i
taft_pKa_plot <- simple_plot_mb_style(taft_pKa_data,
                                      pKa_H2O,
                                      sigma_i,
                                      3,
                                      "black",
                                      expression(paste(
                                        pK[a],
                                        "(",
                                        H[2],
                                        "O)")),
                                      expression(
                                        sigma[I]))

# Regression table for pKa regression
pka_reg_table <- modelsummary(
  lm_pKa_sigma,
  fmt = 3,
  estimate = "{estimate} [{conf.low}, {conf.high}]",
  statistic = NULL,
  gof_omit = "AIC|BIC|Log.Lik.",
  title = "Summary Statistics pKa Model",
                               output = "latex")

# Plot for correlation between sigma i of 
# sulfonyls and alkyls
sulfonyl_alkyl_sigma_plot <- simple_plot_mb_style(
  sigma_sulf,
  sigma_i_alk,
  sigma_i_sulf,
  3,
  "black",
  expression(
    paste(
      sigma[I],
      " of Alkyl Substituents")),
  expression(
    paste(
      sigma[I],
      " of Sulfonyl Substituents")))

# Regression table for pKa regression
sigma_reg_table <- modelsummary(
  lm_sulf_sigma_i,
  fmt = 3,
  estimate = "{estimate} [{conf.low}, {conf.high}]",
  statistic = NULL,
  gof_omit = "AIC|BIC|Log.Lik.",
  title = "Summary Statistics Alkyl & Sulfonyl sigma",
  output = "latex")

#================
# Backbone models
#================
# Model summary plot for IR vs sigma model
model_comparison_plot <- plot_summs(
  lm_backbones_IR,
  lm_sum_back,
  e_nu_lm_backbones,
  omit.coefs = NULL,
  plot.distributions = TRUE,
  model.names = c("v(N-H) Model",
                  "σ(I) Model",
                  "E(HOMO)/ω(E) Model"),
  coefs = c("Intercept" = "(Intercept)",
            "σ para" = "sigma_para_scaled",
            "v(N-H)" = "NH_IR_scaled",
            "Σ(σ(Ι)Ν)" = "sum_s_i_a_scaled",
            "E(HOMO)Nu" = "Nu_scaled",
            "ω(E)" = "Electr_scaled"))

# Model comparison table for IR and sigma models
model_comparison_table <- modelsummary(list(
  "IR Model" = lm_backbones_IR,
  "Sigma i Model" = lm_sum_back,
  "Comp. Descr. Model" = e_nu_lm_backbones),
fmt = 3,
estimate = "{estimate} [{conf.low}, {conf.high}]",
statistic = NULL,
gof_omit = "AIC|BIC|Log.Lik.",
title = "Summary Statistics for Backbone Models",
output = "latex")

# Plot showing the IR values for each backbone
backbone_IR_plot <- simple_plot_mb_style(
  t_half_model,
  Backbone,
  N_H_stretch_IR,
  3,
  "black",
  "Backbone",expression(
    paste(
      "N-H IR Stretching Frequency (", 
      {cm}^{-1},
      ")"))) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1)) +
  scale_x_discrete(labels = parse(
    text = t_half_model$Special_label_exp))

# Plot backbone vs values
t_h_pred_plot_sum <- ggplot(data = t_half_pred_df_2,
                            aes(x = log_half_life,
                                y = 
                                  log_t_half_pred_sum)) +
  geom_abline(slope = 1,
              intercept = 0,
              color = "#4E525A",
              linewidth = 1) +
  geom_point(size = 4,
             fill = "#636C9D",
             pch = 21) +
  geom_label_repel(data = filter(t_half_pred_df_2,
                                 Backbone %in% c(
                                   "CF3-NMs-NMe2",
                                   "F4-NClMs-OMe",
                                   "F4-NTf-OMe")),
                   aes(label = Special_label_exp),
                   nudge_x = -1,
                   nudge_y = 0.6,
                   parse = T,
                   size = 4) +
  geom_label_repel(data = filter(t_half_pred_df_2,
                                 !(Backbone %in% c(
                                   "CF3-NMs-NMe2",
                                   "F4-NClMs-OMe",
                                   "F4-NTf-OMe"))),
                   aes(label = Special_label_exp),
                   nudge_x = 2,
                   nudge_y = -0.4,
                   parse = T,
                   size = 4) +
  labs(x = expression(paste(log[10](t[1/2]),
                            " at 25 \u00B0C")),
       y = expression(paste("Predicted ", 
                            log[10](t[1/2]),
                            " at 25 \u00B0C"))) +
  xlim(0,9) +
  ylim(0,9) +
  themething

#========================
# Alternative descriptors
#========================
# Check colinearity of HOMO and sigma and 
# electr and Hammett
back_HOMO_sigma_plot <- ggplot(hl_e_nu_df) +
  geom_point(aes(x = as.numeric(
    E_HOMO_amido),
    y = sum_s_i_a),
    color = "red",
    size = 3) +
  geom_point(aes(x = as.numeric(
    aryl_electr),
  y = s_i_electr_scaled),
  size = 3) +
  labs(x = "Electr./Nucl. Descriptor",
       y = "Charton-Taft or Hammett param.") +
  themething

# Predicted vs observed plot for alternative descriptor
# model
pred_obs_plot_nu_e_back <- ggplot(data = hl_e_nu_df %>%
                                    mutate(pred_hl =
                            predict(e_nu_lm_backbones)),
                                  aes(x = log_half_life,
                                      y = pred_hl)) +
  geom_abline(slope = 1,
              intercept = 0,
              color = "#4E525A",
              linewidth = 1) +
  geom_point(size = 4,
             color = "#AF0000") +
  geom_label_repel(aes(label = Special_label_exp),
                   nudge_y = 2,
                   parse = T) +
  labs(y = "Predicted log half life",
       x = "Measured log half life") +
  xlim(0,10) +
  ylim(0,10) +
  themething
