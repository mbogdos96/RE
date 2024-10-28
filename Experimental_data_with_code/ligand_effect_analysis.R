library(dplyr)
library(ggplot2)
library(ggrepel)
library(jtools)
library(AICcmodavg)
library(sjPlot)
library(modelsummary)
library(patchwork)

source("SW_plots.R")
source("common_functions.R")
source("eyring.R")
source("Eox_descriptors.R")

#=============
#=============
# Calculations
#=============
#=============

#=========
# Barriers
#=========
# Import the DFT energy data
dft_energies <- read.csv(file.path(
  "./barriers/structure_energies.csv"))

# Caclculate the barriers from the DFT data
# This df used to create DFT data table in SI
barriers_df <- dft_energies %>%
  dplyr::select(Ligand,
                Backbone,
                qh.G.T._SPC,
                State,
                L_type) %>%
  group_by(Ligand,
           Backbone) %>%
  mutate(Barrier_kcal_mol = 
           abs(diff(qh.G.T._SPC))*627.509,
         Temperature = 298.15) %>%
  ungroup()

# Translate into logk_rt to use in modeling
barriers_L <- barriers_df  %>%
  filter(State == "TS",
         Backbone == "F4-NClMs-OMe") %>%
  dplyr::select(Ligand,
                Barrier_kcal_mol,
                Temperature,
                L_type) %>%
  mutate(rate_constant = 2.1e10*Temperature*
           exp(-(Barrier_kcal_mol*4184)/
                 (8.314*Temperature)),
         rate_constant_rt = if_else(
           Temperature != 298.15, 
                                 2.1e10*
                                   Temperature*
                      exp((Temperature/298.15)*
                      log(rate_constant/
                         (2.1e10*Temperature))), 
                                 rate_constant),
         logk_rt = log10(rate_constant_rt))

# Import the buried volume data
Boltz_sterics <- read.csv(file.path(
"./ligand_descriptors/steriprop_boltzmanndistributed.csv"
),
                          sep = ";")

# Merge the DFT data df with the oxidation potential 
# data from SWV
barriers_and_E_L <- merge(barriers_L,
                          echem_data) %>%
  merge(Boltz_sterics) %>%
  mutate(Special_label_exp = sapply(
    Ligand, 
    replace_ligand_expression_og))

# EDA data
EDA_df <- read.csv(file.path(
  "./EDA/EDA_rawvalues.csv"))

# Calculate differences in energies
EDA_delta_GS_TS <- EDA_df %>%
  group_by(Ligand) %>%
  mutate(delta_E_Pauli = E_Pauli - lag(E_Pauli),
         delta_E_elst = E_elst - lag(E_elst),
         delta_E_orb = E_orb - lag(E_orb),
         delta_E_disp = E_disp - lag(E_disp)) %>%
  ungroup() %>%
  filter(Structure == "TS") %>%
  dplyr::select(Ligand,
                delta_E_Pauli,
                delta_E_elst,
                delta_E_orb,
                delta_E_disp,
                L_type) %>%
  mutate(delta_E_sterics = 
           delta_E_Pauli + delta_E_elst,
         sum_sterics_disp = 
           delta_E_sterics + delta_E_disp,
         Special_label_exp = sapply(
           Ligand, 
           replace_ligand_expression_og))

# Bring Vbur data into EDA df
EDA_Vbur_df <- merge(Boltz_sterics,
                     EDA_delta_GS_TS)

#=======
#=======
# Models
#=======
#=======
# Scale and center the data which you will need
k_E_L_scaled <- barriers_and_E_L %>%
  mutate(Vbur_scaled = scale(Vbur,
                             center = T,
                             scale = T),
         E_scaled = scale(Potential,
                          center = T,
                          scale = T),
         L_type = as.factor(L_type))

# Model of only electronics
E_only_model <- lm(data = k_E_L_scaled, 
                   logk_rt ~ E_scaled)

# Model of only electronics and sterics
E_Vbur_model <- lm(data = k_E_L_scaled,
                   logk_rt ~ E_scaled + Vbur_scaled)

# Model where the categorical parameter 
# of coord no is included
E_Vbur_geom_model <- lm(data = k_E_L_scaled,
                        logk_rt ~ E_scaled + 
                          Vbur_scaled + L_type)

# Akaike information criterion for models
model_names <- c("E",
                 "E and S",
                 "E, S and CN")

# Create a list which contains the models
rate_models <- list(E_only_model,
                    E_Vbur_model,
                    E_Vbur_geom_model)

models_table_stats <- modelsummary(rate_models,
  title = 
    "Summary statistics for the Ancillary Ligand Models",
  fmt = 3,
  coef_rename = c("Intercept",
                  "E(ox)",
                  "%Vbur",
                  "Coord. No.") ,
  estimate = "{estimate} [{conf.low}, {conf.high}]",
  statistic = NULL,
  gof_omit = "AIC|BIC|Log.Lik.",
  output = "latex")

# Table comparing models w/AIC
aic_compare_models <- aictab(cand.set = rate_models,
                             modnames = model_names)

# Table to comparing models w/BIC
BIC_compare_models <- data.frame(
  "Model" = c("1 (E)",
              "2 (E & S)",
              "3 (E, S & CN)"),
  "BIC" = c(BIC(E_only_model),
            BIC(E_Vbur_model),
            BIC(E_Vbur_geom_model)))

# Anova for models
# Models 1 and 2
anova_Vbur <- anova(E_only_model,
                    E_Vbur_model)

# Model 3
anova_geom <- anova(E_Vbur_model,
                    E_Vbur_geom_model)

# EDA models
lm_EDA_all <- lm(data = EDA_Vbur_df %>%
                   filter(Ligand != "Fdppe"),
                 sum_sterics_disp ~ Vbur)

lm_EDA_seg <- lm(data = EDA_Vbur_df %>%
                   filter(Ligand != "Fdppe"),
                 sum_sterics_disp ~ Vbur + L_type)

lm_EDA_mono <- lm(data = EDA_Vbur_df %>%
                    filter(Ligand != "Fdppe"
                           & L_type == "monodentate"),
                  sum_sterics_disp ~ Vbur)

lm_EDA_bi <- lm(data = EDA_Vbur_df %>%
                    filter(Ligand != "Fdppe"
                           & L_type == "bidentate"),
                  sum_sterics_disp ~ Vbur)

EDA_pred_df <- EDA_Vbur_df %>%
  mutate(pred_sum = predict(lm_EDA_all,
                            EDA_Vbur_df))
#============
# Predictions
#============
# Create the dataframe with the predictions in the model
logk_pred_df <- k_E_L_scaled %>%
  mutate(logk_pred_no_geom = predict(E_Vbur_model,
                                     k_E_L_scaled),
         logk_pred_geom = predict(E_Vbur_geom_model,
                                  k_E_L_scaled),
         logk_pred_e = predict(E_only_model,
                               k_E_L_scaled),
         Special_label_exp = sapply(
           Ligand, 
           replace_ligand_expression_og))

#======
#======
# Plots
#======
#======

#=========
# Barriers
#=========
# Plot of log rate constant at room temperature 
# vs first oxidation potential
logk_E_plot <- ggplot(data = barriers_and_E_L,
                      aes(x = logk_rt,
                          y = Potential)) +
  geom_point(aes(color = L_type),
             size = 3) +
  geom_label_repel(aes(label = Special_label_exp),
                   nudge_x = 1,
                   nudge_y = -0.02,
                   parse = TRUE) +
  labs(x = "log(k) at 25 \u00B0C",
       y = expression(paste(
         "First Oxidation Potential vs. Fc/", 
                            {"Fc"}^{"+"}, 
                            " (V)"))) +
  themething +
  theme(legend.position = "none")

# EDA Plot
EDA_plot <- ggplot(data = EDA_Vbur_df,
                   aes(x = Vbur,
                       y = sum_sterics_disp)) +
  geom_label_repel(aes(label = Special_label_exp),
                   parse = T,
                   nudge_x = 1) +
  geom_point(data = . %>%
               filter(L_type == "monodentate"),
             fill = "#636C9D",
             size = 4,
             pch = 21) +
  geom_point(data = . %>%
             filter(L_type == "bidentate"),
             fill = "#B71B1B",
           size = 4,
           pch = 21) +
  labs(x = expression(paste("%",
                            V[bur])),
       y = expression(paste(Delta,
                            E[sterics],
                            " + ",
                            Delta,
                            E[disp.],
                            " (kcal/mol)"))) +
    themething

# Plot for comparing the models of electronics, 
# sterics, geometry
models_rates_comparison <- plot_summs(
  E_only_model,
  E_Vbur_model,
  E_Vbur_geom_model,
  omit.coefs = NULL,
  plot.distributions = TRUE,
  model.names = c("Model 1",
                  "Model 2",
                  "Model 3"),
  coefs = c(
    "Intercept" = 
      "(Intercept)",
    "E(ox)" = 
    
      "E_scaled",
    "%Vbur" = 
      "Vbur_scaled",
    "Coord. No." = 
      "L_typemonodentate"))

# Make the plot with all models on it
# Function for making the plots
mk_anc_L_pred_plot <- function(inp_df,
                               nudge_x_1,
                               nudge_y_1,
                               nudge_x_2,
                               nudge_y_2,
                               pred_logk){
  # Create the plot
  out_plot <- ggplot(data = inp_df) +
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
    geom_label_repel(data = filter(
      inp_df,
      (logk_rt/{{pred_logk}}) >= 1),
      aes(x = logk_rt,
          y = {{pred_logk}},
          label = Special_label_exp),
      nudge_x = nudge_x_1,
      nudge_y = nudge_y_1,
      parse = T,
      size = 4) +
    geom_label_repel(data = filter(
      inp_df,
      (logk_rt/{{pred_logk}}) < 1),
      aes(x = logk_rt,
          y = {{pred_logk}},
          label = Special_label_exp),
      nudge_x = nudge_x_2,
      nudge_y = nudge_y_2,
      parse = T, 
      size = 4) +
    geom_point(data = filter(inp_df,
                             L_type == "bidentate"),
               aes(x = logk_rt,
                   y = {{pred_logk}}),
               size = 4,
               fill = "#B71B1B",
               pch = 21) +
    geom_point(data = filter(inp_df,
                             L_type == "monodentate"),
               aes(x = logk_rt,
                   y = {{pred_logk}}),
               size = 5,
               fill = "#636C9D",
               pch = 21) +
    labs(
      x = expression(paste(log[10](k),
                           " at 25 \u00B0C")),
      y = expression(paste("Predicted ", 
                           log[10](k),
                           " at 25 \u00B0C"))) +
    xlim(-20,0) +
    ylim(-20,0) +
    themething
  
  # Output the plot
  return(out_plot)
}

# Plot for model 1
pred_plot_mod1 <- mk_anc_L_pred_plot(logk_pred_df,
                                     -1.3,0.6,1.1,-1,
                                     logk_pred_e)

# Plot for model 2
pred_plot_mod2 <- mk_anc_L_pred_plot(logk_pred_df,
                                     -1.2,0.6,1.1,-0.8,
                                     logk_pred_no_geom)

# Plot for model 3
pred_plot_mod3 <- mk_anc_L_pred_plot(logk_pred_df,
                                     -1.3,0.6,1.1,-1,
                                     logk_pred_geom)

# Combine all three plots
pred_plot_all_logk <- pred_plot_mod1 +
  pred_plot_mod2 +
  pred_plot_mod3 +
  plot_layout(ncol = 3,
              nrow = 1,
              axis_titles = "collect") &
  theme(plot.title = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)))

# Make the plot with only the best model 
# on it for the manuscript
pred_plot_logk_ms <- ggplot() +
  geom_abline(slope = 1,
              intercept = 0,
              color = "#4E525A",
              linewidth = 1) +
  geom_label_repel(data = filter(logk_pred_df,
            (logk_rt/logk_pred_no_geom) >= 1),
                   aes(x = logk_rt,
                       y = logk_pred_no_geom,
                       label = Special_label_exp),
                   nudge_x = -1.3,
                   nudge_y = 0.6,
                   parse = T,
                   size = 4) +
  geom_label_repel(data = filter(logk_pred_df,
                 (logk_rt/logk_pred_no_geom) < 1
                                 & logk_rt > -11),
                   aes(x = logk_rt,
                       y = logk_pred_no_geom,
                       label = Special_label_exp),
                   nudge_x = 1.1,
                   nudge_y = -1,
                   parse = T, 
                   size = 4) +
  geom_label_repel(data = filter(logk_pred_df,
                                 Ligand == "Fdppe"),
                   aes(x = logk_rt,
                       y = logk_pred_no_geom,
                  label = Special_label_exp),
                   nudge_x = -0.5,
                   nudge_y = -2,
                   parse = T, 
                   size = 4) +
  geom_point(data = filter(logk_pred_df,
                           L_type == "bidentate"),
             aes(x = logk_rt,
                 y = logk_pred_no_geom),
             size = 5,
             color = "#B71B1B") +
  geom_point(data = filter(logk_pred_df,
                           L_type == "monodentate"),
             aes(x = logk_rt,
                 y = logk_pred_no_geom),
             size = 5,
             color = "#636C9D") +
  labs(x = expression(paste(log[10](k),
                            " at 25 \u00B0C")),
       y = expression(paste("Predicted ", 
                            log[10](k),
                            " at 25 \u00B0C"))) +
  xlim(-20, 0) +
  ylim(-20,0) +
  themething +
  theme(axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.5)))
