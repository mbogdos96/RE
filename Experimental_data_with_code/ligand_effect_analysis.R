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

#=============
#=============
# Calculations
#=============
#=============

#============
# Descriptors
#============
# Keep only the oxidation potentials and the ligand for simplicity
echem_data <- first_oxidation_potential_all %>%
  dplyr::select(Ligand,
                Potential)

# Import the CEP data
CEP_ligands_ML <- read.csv(file.path("./ligand_descriptors/ML_TEP.csv"))

#Import DFT calculated frequencies
TEP_DFT <- read.csv(file.path("./ligand_descriptors/DFT_TEP.csv"))

# Create a dataframe which contains the DFT and ML TEPs
CO_stretch_method_df <- merge(CEP_ligands_ML,
                              TEP_DFT) %>%
  filter(L_type == "monodentate") %>%
  dplyr::select(Ligand,
         TolmanPrediction,
         Ni_CO_stretch) %>%
  mutate(Ni_CO_stretch = as.numeric(Ni_CO_stretch),
         rel_ML = TolmanPrediction - min(TolmanPrediction),
         rel_DFT = Ni_CO_stretch - min(Ni_CO_stretch),
         Special_label_exp = sapply(Ligand, replace_ligand_expression_og))

# Create dataframe which contains both CO stretch and oxidation potentials
Ni_Rh_TEP_E <- merge(TEP_DFT,
                     first_oxidation_potential_all) %>%
               mutate(Ni_CO_stretch = as.numeric(Ni_CO_stretch),
                      Rh_CO_stretch = as.numeric(Rh_CO_stretch),
                      rel_Ni_CO_stretch = Ni_CO_stretch - min(Ni_CO_stretch,
                                                              na.rm = T),
                      rel_Rh_CO_stretch = Rh_CO_stretch - min(Rh_CO_stretch,
                                                              na.rm = T),
                      Special_label_exp = sapply(Ligand, 
                                                 replace_ligand_expression_og),
                      rel_stretch = ifelse(!is.na(rel_Ni_CO_stretch),
                                           rel_Ni_CO_stretch,
                                           rel_Rh_CO_stretch))

# Check correlation CO stretched and oxidation potential for all, 
# or each category
corr_TEP_E_all <- summary(lm(Potential ~ rel_stretch,
                     data = Ni_Rh_TEP_E))

corr_TEP_E_mono <- summary(lm(Potential ~ Ni_CO_stretch,
                             data = Ni_Rh_TEP_E))

corr_TEP_E_bi <- summary(lm(Potential ~ rel_Rh_CO_stretch,
                             data = Ni_Rh_TEP_E))

# Pd Mulliken charges and HOMO energies
Pd_charges <- read.csv(file.path("./complex_descriptors/HOMO_and_charges.csv")) %>%
  merge(echem_data) %>%
  mutate(HOMO_energy = as.numeric(HOMO_energy),
  Special_label_exp = sapply(Ligand, 
                             replace_ligand_expression_og))

lm_HOMO_E <- lm(Potential ~ HOMO_energy,
                data = Pd_charges)

lm_Mulliken_E <- lm(Potential ~ Pd_Mulliken_charge,
                    data = Pd_charges)

Mull_plot <- ggplot() +
  geom_smooth(data = filter(Pd_charges,
                            L_type == "bidentate"),
              aes(x = Potential,
                  y = Pd_Mulliken_charge),
              method = lm,
              se = T,
              alpha = 0.3,
              color = "#B71B1B") +
  geom_smooth(data = filter(Pd_charges,
                            L_type == "monodentate"),
              aes(x = Potential,
                  y = Pd_Mulliken_charge),
              method = lm,
              se = T,
              alpha = 0.3,
              color = "#636C9D") +
  geom_label_repel(data = filter(Pd_charges,
                                 L_type == "bidentate"),
                   aes(x = Potential,
                       y = Pd_Mulliken_charge,
                       label = Special_label_exp),
                   parse = T,
                   nudge_y = -0.1,
                   nudge_x = 0.01,
                   size = 3) +
  geom_label_repel(data = filter(Pd_charges,
                                 L_type == "monodentate"),
                   aes(x = Potential,
                       y = Pd_Mulliken_charge,
                       label = Special_label_exp),
                   parse = T,
                   nudge_y = 0.1,
                   nudge_x = -0.01,
                   size = 3) +
  geom_point(data = filter(Pd_charges,
                           L_type == "bidentate"),
             aes(x = Potential,
                 y = Pd_Mulliken_charge),
             color = "#B71B1B",
             size = 4) +
  geom_point(data = filter(Pd_charges,
                           L_type == "monodentate"),
             aes(x = Potential,
                 y = Pd_Mulliken_charge),
             color = "#636C9D",
             size = 4) +
  labs(x = expression(paste("Potential (V, Fc/", 
                            {Fc^{"+"}} ,
                            ")")),
       y = "Mulliken Charge at Pd") +
  themething

HOMO_plot <- ggplot() +
  geom_smooth(data = filter(Pd_charges,
                            L_type == "bidentate"),
              aes(x = Potential,
                  y = HOMO_energy),
              method = lm,
              se = T,
              alpha = 0.3,
              color = "#B71B1B") +
  geom_smooth(data = filter(Pd_charges,
                            L_type == "monodentate"),
              aes(x = Potential,
                  y = HOMO_energy),
              method = lm,
              se = T,
              alpha = 0.3,
              color = "#636C9D") +
  geom_label_repel(data = filter(Pd_charges,
                                 L_type == "bidentate"),
                   aes(x = Potential,
                       y = HOMO_energy,
                       label = Special_label_exp),
                   parse = T,
                   nudge_y = -0.002,
                   nudge_x = -0.01,
                   size = 3) +
  geom_label_repel(data = filter(Pd_charges,
                                 L_type == "monodentate"),
                   aes(x = Potential,
                       y = HOMO_energy,
                       label = Special_label_exp),
                   parse = T,
                   nudge_y = 0.002,
                   nudge_x = 0.01,
                   size = 3) +
  geom_point(data = filter(Pd_charges,
                           L_type == "bidentate"),
             aes(x = Potential,
                 y = HOMO_energy),
             color = "#B71B1B",
             size = 4) +
  geom_point(data = filter(Pd_charges,
                           L_type == "monodentate"),
             aes(x = Potential,
                 y = HOMO_energy),
             color = "#636C9D",
             size = 4) +
  labs(x = expression(paste("Potential (V, Fc/", 
                            {Fc^{"+"}} ,
                            ")")),
       y = "HOMO Energy") +
  themething

desc_E_corr_plot <- Mull_plot +
  HOMO_plot + 
  plot_layout(guides = "collect",
              axes = "collect") &
  theme(axis.title = element_text(size = rel(1.2)),
        axis.text = element_text(size = rel(1.2)))

#=========
# Barriers
#=========
# Import the DFT energy data
dft_energies <- read.csv(file.path("./barriers/structure_energies.csv"))

barriers_df <- dft_energies %>%
  dplyr::select(Ligand,
                Backbone,
                qh.G.T._SPC,
                State,
                L_type) %>%
  group_by(Ligand,
           Backbone) %>%
  mutate(Barrier_kcal_mol = abs(diff(qh.G.T._SPC))*627.509,
         Temperature = 298.15) %>%
  ungroup()
  

barriers_L <- barriers_df  %>%
  filter(State == "TS",
         Backbone == "F4-NClMs-OMe") %>%
  dplyr::select(Ligand,
                Barrier_kcal_mol,
                Temperature,
                L_type) %>%
  mutate(rate_constant = 2.1e10*Temperature*exp(-(Barrier_kcal_mol*4184)/
                                                  (8.314*Temperature)),
         rate_constant_rt = if_else(Temperature != 298.15, 
                                 2.1e10*
                                   Temperature*
                                   exp((Temperature/298.15)*
                                         log(rate_constant/
                                               (2.1e10*Temperature))), 
                                 rate_constant),
         logk_rt = log10(rate_constant_rt))

# Import the buried volume data
Boltz_sterics <- read.csv(file.path(
  "./ligand_descriptors/steriprop_boltzmanndistributed.csv"),
                          sep = ";")

# Merge the DFT data df with the oxidation potential data from SWV
barriers_and_E_L <- merge(barriers_L,
                          echem_data) %>%
  merge(Boltz_sterics) %>%
  mutate(Special_label_exp = sapply(Ligand, 
                                    replace_ligand_expression_og))

b_E_L_add_HOMO <- merge(barriers_and_E_L,
                        Pd_charges) %>%
  mutate(Vbur_scaled = scale(Vbur,
                             center = T,
                             scale = T),
         HOMO_scaled = scale(HOMO_energy,
                             scale = T,
                             center = T))

lm_HOMO_logk <- lm(data = b_E_L_add_HOMO,
                   logk_rt ~ HOMO_scaled + Vbur_scaled)

# EDA data
EDA_df <- read.csv(file.path(
  "./EDA/EDA_rawvalues.csv"))

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
  mutate(delta_E_sterics = delta_E_Pauli + delta_E_elst,
         sum_sterics_disp = delta_E_sterics + delta_E_disp,
         Special_label_exp = sapply(Ligand, 
                                    replace_ligand_expression_og))
  
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
E_only_model <- lm(data = k_E_L_scaled, logk_rt ~ E_scaled)

# Model of only electronics and sterics
E_Vbur_model <- lm(data = k_E_L_scaled,
                   logk_rt ~ E_scaled + Vbur_scaled)

# Model where the categorical parameter of coord no is included
E_Vbur_geom_model <- lm(data = k_E_L_scaled,
                        logk_rt ~ E_scaled + Vbur_scaled + L_type)

# Akaike information criterion for models
model_names <- c("E",
                 "E and S",
                 "E, S and CN")

rate_models <- list(E_only_model,
                    E_Vbur_model,
                    E_Vbur_geom_model)

models_table_stats <- modelsummary(rate_models,
  title = "Summary statistics for the 3 models unders study",
  fmt = 3,
  coef_rename = c("Intercept",
                  "Oxidation Potential",
                  "Buried Volume",
                  "Coordination Number") ,
  estimate = "{estimate} [{conf.low}, {conf.high}]",
  statistic = NULL,
  gof_omit = "AIC|BIC|Log.Lik.",
  output = "latex")

aic_compare_models <- aictab(cand.set = rate_models,
                             modnames = model_names)

BIC_compare_models <- data.frame(
  "Model" = c("E",
              "E and S",
              "E, S and CN"),
  "BIC" = c(BIC(E_only_model),
            BIC(E_Vbur_model),
            BIC(E_Vbur_geom_model)))

# Anova for models
anova_Vbur <- anova(E_only_model,
                    E_Vbur_model)

anova_geom <- anova(E_Vbur_model,
                    E_Vbur_geom_model)

# Check what the outcome is for mixing experimental with DFT barriers
eyring_rate_const_df <- eyring_summary %>%
  mutate(rate_constant_rt = if_else(Temperature_K != 298.15, 
                                    2.1e10*
                                      Temperature_K*
                                      exp((Temperature_K/298.15)*
                                            log(rate_constant/
                                                  (2.1e10*Temperature_K))), 
                                    rate_constant))

mixed_bar <- barriers_and_E_L %>%
  mutate(exp_rate_constant = 
           c(NA,NA,NA,1.301464e-08,
             9.451181e-05,NA,4.281247e-07,3.882400e-05,7.248984e-05),
         mixed_logkrt = if_else(is.na(exp_rate_constant ),
                                logk_rt,
                                exp_rate_constant))

mixed_E <- lm(data = mixed_bar,
                mixed_logkrt ~ Potential)

mixed_E_S <- lm(data = mixed_bar,
                mixed_logkrt ~ Potential + Vbur)

mixed_E_S_CN <- lm(data = mixed_bar,
                mixed_logkrt ~ Potential + Vbur + L_type)

mixed_rate_models <- list(mixed_E,
                    mixed_E_S,
                    mixed_E_S_CN)

aic_compare_mixed_models <- aictab(cand.set = mixed_rate_models,
                             modnames = model_names)

BIC_compare_mixed_models <- list(c(BIC(mixed_E),
                                   BIC(mixed_E_S),
                                   BIC(mixed_E_S_CN)))

anova_mixed_S <- anova(mixed_E,
                       mixed_E_S)

anova_mixed_CN <- anova(mixed_E_S,
                        mixed_E_S_CN)

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

#===================
# Mediation analysis
#===================
# Checking relative magnitudes of direct and indirect effects
lm_S_indirect_scaled <- lm(E_scaled ~ L_type 
                           + Vbur_scaled,
                           data = k_E_L_scaled)

lm_E_S_direct_raw <- lm(logk_rt ~ Potential + Vbur,
                        data = k_E_L_scaled)

lm_S_indirect_raw <- lm(Potential ~ L_type 
                        + Vbur,
                        data = k_E_L_scaled)

rel_S_effect_df <- data.frame("Scaled Values" = c(T,T,T,
                                                  F,F,F),
                              "Effect Size" = c(
                                E_Vbur_model$coefficients[[3]],
                                lm_S_indirect_scaled$coefficients[[3]]
                                *E_Vbur_model$coefficients[[2]],
                                lm_S_indirect_scaled$coefficients[[2]]
                                *E_Vbur_model$coefficients[[2]],
                                lm_E_S_direct_raw$coefficients[[3]],
                                lm_S_indirect_raw$coefficients[[3]]
                                *lm_E_S_direct_raw$coefficients[[2]],
                                lm_S_indirect_raw$coefficients[[2]]
                                *lm_E_S_direct_raw$coefficients[[2]]
                              ),
                              "Effect Type" = c("direct",
                                               "through electronics",
                                               "through CN",
                                               "direct",
                                               "through electronics",
                                               "through CN"))

total_S_effect <- data.frame("Scaled Values" = c(T,F),
                             "S_Total_Effect" = c(
                               sum((rel_S_effect_df %>%
                                      filter(Scaled.Values == T))$Effect.Size),
                               sum((rel_S_effect_df %>%
                                      filter(Scaled.Values == F))$Effect.Size)
                             )) %>%
  mutate(Electronic_Effect = c(E_Vbur_model$coefficients[[2]],
                                 lm_E_S_direct_raw$coefficients[[2]]),
         "Electronic_v_Steric_Fold_D" = Electronic_Effect/S_Total_Effect)

mediation_mod_list <- list(E_Vbur_model,
                           lm_S_indirect_scaled,
                           lm_S_indirect_raw,
                           lm_E_S_direct_raw)

mediation_table_stats <- modelsummary(
  mediation_mod_list,
  title = "Model statistics for models used in mediation analysis",
  fmt = 3,
  estimate = "{estimate} [{conf.low}, {conf.high}]",
  statistic = NULL,
  gof_omit = "AIC|BIC|Log.Lik.",
  output = "latex")

# Try calculating effects from lm w/std errors to make error bars
med_effects_errors_df <- summary(E_Vbur_model)$coefficients %>%
  rbind(summary(lm_S_indirect_scaled)$coefficients) %>%
  as.data.frame() %>%
  rename(Std_error = "Std. Error") %>%
  mutate(effect = c("none",
                    "E",
                    "S_D",
                    "none",
                    "CN_E",
                    "S_E"),
         upper_est = Estimate + Std_error,
         lwr_est = Estimate - Std_error) %>%
  dplyr::select(Estimate, effect, upper_est, lwr_est) %>%
  mutate(scaling = T)

med_effects_errors_raw_df <- summary(lm_E_S_direct_raw)$coefficients %>%
  rbind(summary(lm_S_indirect_raw)$coefficients) %>%
  as.data.frame() %>%
  rename(Std_error = "Std. Error") %>%
  mutate(effect = c("none",
                    "E",
                    "S_D",
                    "none",
                    "CN_E",
                    "S_E"),
         upper_est = Estimate + Std_error,
         lwr_est = Estimate - Std_error) %>%
  dplyr::select(Estimate, effect, upper_est, lwr_est) %>%
  mutate(scaling = F)

mediation_estimates_errors_df <- data.frame(
  "Estimate" = c(med_effects_errors_df$Estimate[[2]],
                 med_effects_errors_df$Estimate[[3]],
                 med_effects_errors_df$Estimate[[2]]
                 *med_effects_errors_df$Estimate[[6]],
                 med_effects_errors_df$Estimate[[2]]
                 *med_effects_errors_df$Estimate[[5]],
                 med_effects_errors_df$Estimate[[3]] +
                 med_effects_errors_df$Estimate[[2]]
                 *med_effects_errors_df$Estimate[[6]] +
                 med_effects_errors_df$Estimate[[2]]
                 *med_effects_errors_df$Estimate[[5]]),
  "Upper" = c(med_effects_errors_df$upper_est[[2]],
              med_effects_errors_df$upper_est[[3]],
              med_effects_errors_df$upper_est[[2]]
              *med_effects_errors_df$upper_est[[6]],
              med_effects_errors_df$upper_est[[2]]
              *med_effects_errors_df$upper_est[[5]],
              med_effects_errors_df$upper_est[[3]] +
              med_effects_errors_df$upper_est[[2]]
              *med_effects_errors_df$upper_est[[6]] +
              med_effects_errors_df$upper_est[[2]]
              *med_effects_errors_df$upper_est[[5]]),
  "Lower" = c(med_effects_errors_df$lwr_est[[2]],
              med_effects_errors_df$lwr_est[[3]],
              med_effects_errors_df$lwr_est[[2]]
              *med_effects_errors_df$lwr_est[[6]],
              med_effects_errors_df$lwr_est[[2]]
              *med_effects_errors_df$lwr_est[[5]],
              med_effects_errors_df$lwr_est[[3]] + 
                med_effects_errors_df$lwr_est[[2]]
              *med_effects_errors_df$lwr_est[[6]] + 
                med_effects_errors_df$lwr_est[[2]]
              *med_effects_errors_df$lwr_est[[5]]),
  "Effect" = c("E",
               "S_D",
               "S_I",
               "S_CN",
               "S_tot"))

mediation_estimates_errors_raw_df <- data.frame(
  "Estimate" = c(med_effects_errors_raw_df$Estimate[[2]],
                 med_effects_errors_raw_df$Estimate[[3]],
                 med_effects_errors_raw_df$Estimate[[2]]*
                   med_effects_errors_raw_df$Estimate[[6]],
                 med_effects_errors_raw_df$Estimate[[2]]*
                   med_effects_errors_raw_df$Estimate[[5]],
                 med_effects_errors_raw_df$Estimate[[3]] +
                 med_effects_errors_raw_df$Estimate[[2]]*
                   med_effects_errors_raw_df$Estimate[[6]] +
                 med_effects_errors_raw_df$Estimate[[2]]*
                   med_effects_errors_raw_df$Estimate[[5]]),
  "Upper" = c(med_effects_errors_raw_df$upper_est[[2]],
              med_effects_errors_raw_df$upper_est[[3]],
              med_effects_errors_raw_df$upper_est[[2]]*
                med_effects_errors_raw_df$upper_est[[6]],
              med_effects_errors_raw_df$upper_est[[2]]*
                med_effects_errors_raw_df$upper_est[[5]],
              med_effects_errors_raw_df$upper_est[[3]] +
              med_effects_errors_raw_df$upper_est[[2]]*
                med_effects_errors_raw_df$upper_est[[6]] +
              med_effects_errors_raw_df$upper_est[[2]]*
                med_effects_errors_raw_df$upper_est[[5]]),
  "Lower" = c(med_effects_errors_raw_df$lwr_est[[2]],
              med_effects_errors_raw_df$lwr_est[[3]],
              med_effects_errors_raw_df$lwr_est[[2]]*
                med_effects_errors_raw_df$lwr_est[[6]],
              med_effects_errors_raw_df$lwr_est[[2]]*
                med_effects_errors_raw_df$lwr_est[[5]],
              med_effects_errors_raw_df$lwr_est[[3]] +
              med_effects_errors_raw_df$lwr_est[[2]]*
                med_effects_errors_raw_df$lwr_est[[6]] +
              med_effects_errors_raw_df$lwr_est[[2]]*
                med_effects_errors_raw_df$lwr_est[[5]]),
  "Effect" = c("E",
               "S_D",
               "S_I",
               "S_CN",
               "S_tot"))

med_ef_df_all <- mediation_estimates_errors_raw_df %>%
  rbind(mediation_estimates_errors_df) %>%
  mutate(scaling = c(F,F,F,F,F,T,T,T,T,T))

mediation_estimates_plot <- ggplot(data = med_ef_df_all) +
  geom_vline(xintercept = 0, 
              linetype = "dashed", 
              color = "#4E525A", 
              linewidth = 1) +
  geom_errorbarh(aes(xmax = Upper,
                     xmin = Lower,
                     y = Effect,
                     color = scaling),
                 linewidth = 0.75,
                 alpha = 0.3) +
  geom_point(aes(x = Estimate,
                 y = Effect,
                 color = scaling),
             size = 5) +
  scale_color_manual(values = c("#636C9D","#B71B1B")) +
  scale_y_discrete(labels = parse(text = 
                                    c("e",
                                      "s^{cn}",
                                      "s^{d}",
                                      "{s}^{i}",
                                      "s^{tot}"))) +
  themething +
  theme(legend.position = "none")

# Exploratory only - weird finding
lm_interaction <- lm(logk_rt ~ E_scaled + 
                       Vbur_scaled +
                       E_scaled:Vbur_scaled,
                     data = k_E_L_scaled)

lm_med_as_int <- lm(logk_rt ~ E_scaled + 
                      Vbur_scaled +
                      E_scaled:Vbur_scaled +
                      E_scaled:L_type,
                    data = k_E_L_scaled)

pred_df <- k_E_L_scaled %>%
  mutate(pred_rate = predict(lm_interaction,
                             k_E_L_scaled),
         pred_rate_full = predict(lm_med_as_int,
                                  k_E_L_scaled))

pred_plot_el <- ggplot(data = pred_df) +
  geom_abline(slope = 1,
              intercept = 0,
              color = "#4E525A",
              linewidth = 1) +
  geom_point(aes(x = logk_rt,
                 y = pred_rate)) +
  geom_point(aes(x = logk_rt,
                 y = pred_rate_full),
             color = "red")

#======
#======
# Plots
#======
#======

#============
# Descriptors
#============
# Plot comparing DFT vs ML TEP
TEP_DFT_ML_plot <- ggplot(data = CO_stretch_method_df) +
  geom_label_repel(aes(x = Ni_CO_stretch,
                       y = TolmanPrediction,
                       label = Special_label_exp),
                   parse = T,
                   nudge_x = 0.5,
                   nudge_y = 1,
                   size = 5) +
  labs(x = expression(paste("DFT CO Stretching Frequency (", 
                            {cm}^{-1}, 
                            ")")),
       y = expression(paste("Machine Learning CO Stretching Frequency (", 
                            {cm}^{-1}, 
                            ")"))) +
  geom_point(aes(x = Ni_CO_stretch,
                 y = TolmanPrediction),
             size = 5) +
  themething

# Plot comparing DFT vs ML TEP relative values
TEP_DFT_ML_plot_rel <- ggplot(data = CO_stretch_method_df) +
  geom_label_repel(aes(x = rel_DFT,
                       y = rel_ML,
                       label = Special_label_exp),
                   nudge_x = 0.05,
                   nudge_y = 0.1,
                   size = 5) +
  labs(x = expression(paste("Relative DFT CO Stretching Frequency (", 
                            {cm}^{-1}, ")")),
       y = expression(paste(
         "Relative Machine Learning CO Stretching Frequency (", 
         {cm}^{-1}, ")"))) +
  geom_point(aes(x = rel_DFT,
                 y = rel_ML),
             size = 5) +
  themething

# Plot for DFT TEP Ni or Rh vs E ox
E_v_TEP_Ni_Rh <- ggplot() +
  geom_ribbon(data = data.frame("Potential" = c(0.65, 
                                                1.4),
                                "Upper" = c(18, 
                                            37),
                                "Lower" = c(0, 
                                            19)),
              aes(x = Potential,
                  ymax = Upper,
                  ymin = Lower),
              alpha = 0.1) +
  annotate(geom = "text",
           x = 0.85,
           y = 30,
           label = "Expected Relationship",
           alpha = 0.3,
           size = 6) +
  geom_label_repel(data = filter(Ni_Rh_TEP_E,
                                 !is.na(rel_Ni_CO_stretch)),
                   aes(x = Potential,
                       y = rel_Ni_CO_stretch,
                       label = Special_label_exp),
                   nudge_x = 0.01,
                   nudge_y = 3,
                   parse = TRUE,
                   size = 5) +
  geom_label_repel(data = filter(Ni_Rh_TEP_E,
                                 !is.na(Rh_CO_stretch)),
                   aes(x = Potential,
                       y = rel_Rh_CO_stretch,
                       label = Special_label_exp),
                   nudge_x = 0.05,
                   nudge_y = -2,
                   parse = TRUE,
                   size = 5) +
  geom_point(data = filter(Ni_Rh_TEP_E,
                           !is.na(rel_Rh_CO_stretch)),
             aes(x = Potential,
                 y = rel_Rh_CO_stretch),
             color = "#B71B1B",
             size = 5) +
  geom_point(data = filter(Ni_Rh_TEP_E,
                           !is.na(rel_Ni_CO_stretch)),
             aes(x = Potential,
                 y = rel_Ni_CO_stretch),
             color = "#636C9D",
             size = 5) +
  labs(x = expression(paste("Palladacycle ", E^{ox} ," vs. Fc/", 
                            {"Fc"}^{"+"}, " (V)")),
       y = expression(paste("Rh(I) ", {v}[CO], " (rel., ",
                            "DFT, ", {cm}^{-1}, ")"))) +
  scale_y_continuous(sec.axis = sec_axis(~., 
                                         name = expression(
                                           paste(
                                             "Ni(0) ", {v}[CO], " (rel.",
                                             ", DFT, ", {cm}^{-1}, ")")))) +
  coord_cartesian(xlim = c(0.7,1.3)) +
  themething +
  theme(axis.title.y.right = element_text(color = "#636C9D"),
        axis.title.y.left = element_text(color = "#B71B1B"),
        axis.title = element_text(size = rel(1.75)),
        axis.text = element_text(size = rel(1.75)))

#=========
# Barriers
#=========
# Plot of log rate constant at room temperature vs first oxidation potential
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
       y = expression(paste("First Oxidation Potential vs. Fc/", 
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
             color = "#B71B1B",
             size = 4) +
  geom_point(data = . %>%
             filter(L_type == "bidentate"),
             color = "#636C9D",
           size = 4) +
  labs(x = expression(paste("%",
                            V[bur])),
       y = expression(paste(Delta,
                            E[sterics],
                            " + ",
                            Delta,
                            E[disp.],
                            " (kcal/mol)"))) +
    themething

#============
# Predictions
#============
# Create the dataframe with the predictions in the model
logk_pred_df <- k_E_L_scaled %>%
  mutate(logk_pred_no_geom = predict(E_Vbur_model,
                                     k_E_L_scaled),
         pred_no_geom_upr = predict(E_Vbur_model,
                                    k_E_L_scaled,
                                    interval = "prediction",
                                    level = 0.33)[,'upr'], 
         pred_no_geom_lwr = predict(E_Vbur_model,
                                    k_E_L_scaled,
                                    interval = "prediction",
                                    level = 0.33)[,'lwr'],
         pred_no_geom_range = pred_no_geom_upr - pred_no_geom_lwr,
         logk_pred_geom = predict(E_Vbur_geom_model,
                                  k_E_L_scaled),
         pred_geom_upr = predict(E_Vbur_geom_model,
                                 k_E_L_scaled,
                                 interval = "prediction",
                                 level = 0.33)[,'upr'], 
         pred_geom_lwr = predict(E_Vbur_geom_model,
                                 k_E_L_scaled,
                                 interval = "prediction",
                                 level = 0.33)[,'lwr'],
         pred_geom_range = pred_geom_upr - pred_geom_lwr,
         logk_pred_e = predict(E_only_model,
                               k_E_L_scaled),
         pred_e_upr = predict(E_only_model,
                              k_E_L_scaled,
                              interval = "prediction",
                              level = 0.33)[,'upr'], 
         pred_e_lwr = predict(E_only_model,
                              k_E_L_scaled,
                              interval = "prediction",
                              level = 0.33)[,'lwr'],
         pred_e_range = pred_e_upr - pred_e_lwr,
         Special_label_exp = sapply(Ligand, 
                                    replace_ligand_expression_og))

# Plot for comparing the models of electronics, sterics, geometry
models_rates_comparison <- plot_summs(E_only_model,
                                      E_Vbur_model,
                                      E_Vbur_geom_model,
                                      omit.coefs = NULL,
                                      plot.distributions = TRUE,
                                      model.names = c("EGS",
                                                      "EGS, SGS, STS",
                                                      "EGS, SGS, STS, OTS"),
                                      coefs = c(
                                        "Intercept" = 
                                          "(Intercept)",
                                        "Oxidation Potential" = 
                                          "E_scaled",
                                        "Buried Volume" = 
                                          "Vbur_scaled",
                                        "Coordination Number" = 
                                          "L_typemonodentate")) +
  theme(legend.position = "none")

# Make the plot with all models on it
pred_plot_mod1 <- ggplot(data = logk_pred_df) +
  geom_abline(slope = 1,
              intercept = 0,
              color = "#4E525A",
              linewidth = 1) +
  geom_label_repel(data = filter(logk_pred_df,
                                 (logk_rt/logk_pred_e) >= 1),
                   aes(x = logk_rt,
                       y = logk_pred_e,
                       label = Special_label_exp),
                   nudge_x = -1.3,
                   nudge_y = 0.6,
                   parse = T,
                   size = 2) +
  geom_label_repel(data = filter(logk_pred_df,
                                 (logk_rt/logk_pred_e) < 1),
                   aes(x = logk_rt,
                       y = logk_pred_e,
                       label = Special_label_exp),
                   nudge_x = 1.1,
                   nudge_y = -1,
                   parse = T, 
                   size = 2) +
  geom_point(data = filter(logk_pred_df,
                           L_type == "bidentate"),
             aes(x = logk_rt,
                 y = logk_pred_e),
             size = 5,
             color = "#B71B1B") +
  geom_point(data = filter(logk_pred_df,
                           L_type == "monodentate"),
             aes(x = logk_rt,
                 y = logk_pred_e),
             size = 5,
             color = "#636C9D") +
  labs(title = "E Model",
    x = expression(paste(log[10](k)," at 25 \u00B0C")),
       y = expression(paste("Predicted ", log[10](k)," at 25 \u00B0C"))) +
  xlim(-20, 0) +
  ylim(-20,0) +
  themething
  
pred_plot_mod2 <- ggplot(data = logk_pred_df) +
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
                   size = 2) +
  geom_label_repel(data = filter(logk_pred_df,
                                 (logk_rt/logk_pred_no_geom) < 1),
                   aes(x = logk_rt,
                       y = logk_pred_no_geom,
                       label = Special_label_exp),
                   nudge_x = 1.1,
                   nudge_y = -1,
                   parse = T, 
                   size = 2) +
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
  labs(title = "E and S Model",
    x = expression(paste(log[10](k)," at 25 \u00B0C")),
       y = expression(paste("Predicted ", log[10](k)," at 25 \u00B0C"))) +
  xlim(-20, 0) +
  ylim(-20,0) +
  themething +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())
  
pred_plot_mod3 <- ggplot(data = logk_pred_df) +
  geom_abline(slope = 1,
              intercept = 0,
              color = "#4E525A",
              linewidth = 1) +
  geom_label_repel(data = filter(logk_pred_df,
                                 (logk_rt/logk_pred_geom) >= 1),
                   aes(x = logk_rt,
                       y = logk_pred_geom,
                       label = Special_label_exp),
                   nudge_x = -1.3,
                   nudge_y = 0.6,
                   parse = T,
                   size = 2) +
  geom_label_repel(data = filter(logk_pred_df,
                                 (logk_rt/logk_pred_geom) < 1),
                   aes(x = logk_rt,
                       y = logk_pred_geom,
                       label = Special_label_exp),
                   nudge_x = 1.1,
                   nudge_y = -1,
                   parse = T, 
                   size = 2) +
  geom_point(data = filter(logk_pred_df,
                           L_type == "bidentate"),
             aes(x = logk_rt,
                 y = logk_pred_geom),
             size = 5,
             color = "#B71B1B") +
  geom_point(data = filter(logk_pred_df,
                           L_type == "monodentate"),
             aes(x = logk_rt,
                 y = logk_pred_geom),
             size = 5,
             color = "#636C9D") +
  labs(title = "E, S and CN Model" ,
    x = expression(paste(log[10](k)," at 25 \u00B0C")),
       y = expression(paste("Predicted ", log[10](k)," at 25 \u00B0C"))) +
  xlim(-20, 0) +
  ylim(-20,0) +
  themething +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

pred_plot_all_logk <- pred_plot_mod1 +
  pred_plot_mod2 +
  pred_plot_mod3 +
  plot_layout(ncol = 3,
              nrow = 1,
              axis_titles = "collect") &
  theme(plot.title = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)))

# Make the plot with only the best model on it for the manuscript
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
  labs(x = expression(paste(log[10](k)," at 25 \u00B0C")),
       y = expression(paste("Predicted ", log[10](k)," at 25 \u00B0C"))) +
  xlim(-20, 0) +
  ylim(-20,0) +
  themething +
  theme(axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.5)))

# Make the plot with both models on it
pred_plot_all_logk_errors <- ggplot(data = logk_pred_df) +
  geom_point(aes(x = factor(Ligand,
                            levels = unique(Ligand)),
                 y = pred_e_range),
             color = "pink",
             size = 3) +
  geom_point(aes(x = factor(Ligand,
                            levels = unique(Ligand)),
                 y = pred_no_geom_range),
             color = "#636C9D",
             size = 3) +
  geom_point(aes(x = factor(Ligand,
                            levels = unique(Ligand)),
                 y = pred_geom_range),
             color = "#B71B1B",
             size = 3) +
  labs(y = "Width of 67% Prediction Interval") +
  scale_x_discrete(labels = parse(
    text = logk_pred_df$Special_label_exp
  )) +
  themething +
  theme(axis.title.x = element_blank())