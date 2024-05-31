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

# Define the molecular weight of the internal standard used in all cases
internal_standard_molecular_weight <- 246.22

#=============
# Calculations
#=============

# Calculations of parameters
taft_ref_data <- read.csv(file.path("./half_life_data/taft_19_f.csv"))

lm_19F_sigma <- lm(data = taft_ref_data,
                   sigma_i ~ delta_F)

# Method fails for CF3, needs to be outlined
F4_aryl <- data.frame(delta_F = c(0.74,-2.63)) 

F4_aryl_pred <- F4_aryl %>%
  mutate(sigma_i = predict(lm_19F_sigma,F4_aryl)) 

# Sigma i for CF3 arene using COOH regression method
taft_pKa_data <- read.csv(file.path("./half_life_data/taft_pKa.csv"))

lm_pKa_sigma <- lm(data = taft_pKa_data,
                   sigma_i ~ pKa_H2O)

bisCF3_aryl <- data.frame(pKa_H2O = c(3.34)) 

bisCF3_aryl_pred <- bisCF3_aryl %>%
  mutate(sigma_i = predict(lm_pKa_sigma,bisCF3_aryl))

# Sigma i for ClMs using interpolation
sigma_sulf <- read.csv(file.path("./half_life_data/taft_sulf_alk.csv"))

lm_sulf_sigma_i <- lm(data = sigma_sulf,
                      sigma_i_sulf ~ sigma_i_alk)

clMs_sigma_i <- data.frame(sigma_i_alk = 0.13) 

clMs_sigma_i_pred <- clMs_sigma_i %>%
  mutate(sigma_i_sulf = predict(lm_sulf_sigma_i,clMs_sigma_i))

# Import the data and then convert values to SI units
half_life_data <- read.csv(
  file.path("./half_life_data/half_lives_backbones_PtBu3.csv")) %>%
                  mutate(Temperature_K = Temperature_C + 273.15) %>%
                  mutate(time_s = Time_min*60) %>%
                  mutate(Internal_standard_mol = 
                           (Internal_standard_mass_mg*0.001)/
                           internal_standard_molecular_weight) %>%
                  mutate(Pd_II_mol= Pd_II_integral*Internal_standard_mol) %>%
                  group_by(Backbone) %>%
                  mutate(Pd_II_mol_norm = Pd_II_mol/max(Pd_II_mol)) %>%
                  ungroup()

# Function which fits first order rate equation to estimate a rate constant
calculate_rate_constant <- function(df,x,y) {
  unique_backbones <- unique(df$Backbone)
  
  result <- data.frame()
  
  for (backbone in unique_backbones) {
    subset_data <- filter(df, Backbone == backbone)
    model_fit <- try(
      nls(Pd_II_mol ~ A * exp(-B * time_s), 
          data = subset_data, 
          start = list(A = x, B = y)
      )
    )
    
    if (!inherits(model_fit, "try-error")) {
      coefficients <- coef(model_fit)
      result <- bind_rows(result, 
                          data.frame(Backbone = backbone, 
                                     A0 = coefficients[["A"]], 
                                     rate_constant = coefficients[["B"]]))
    }
  }
  
  result <- result
  
  return(result)
}

# Apply the function to the backbone data, manipulate df and calculate t half
backbone_half_lives_PtBu3 <- calculate_rate_constant(half_life_data,
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
                                         log(rate_constant/(2.1e10*Temperature_K))), 
                                 rate_constant)
  ) %>%
  mutate(half_life_s = log(2)/rate_constant,
         log_half_life = log10(half_life_s)) %>%
  mutate(sigma_para_scaled = as.numeric(scale(sigma_para,
                                              center = TRUE,
                                              scale = TRUE)),
         sigma_i_n_scaled = as.numeric(scale(sigma_i_n,
                                             center = TRUE,
                                             scale = TRUE)),
         sigma_i_ar_scaled = as.numeric(scale(sigma_i_ar,
                                              center = TRUE,
                                              scale = TRUE))) %>%
  mutate(NH_IR_scaled = scale(as.numeric(N_H_stretch_IR),
                              center = TRUE,
                              scale = TRUE))

# Remove NH from data before building model, as it is outwith measurement 
# range for rt NMR
t_half_model <- backbone_half_lives_PtBu3 %>%
  mutate(logk = log10(rate_constant),
         Backbone = na_if(Backbone,"CF3-NH-Ph"),
         Special_label_exp = sapply(Backbone,
                                    replace_backbone_labels)) %>%
  na.omit()

# Create a model that relates the IR and sigma para to the log of the half life
lm_backbones_IR <- lm(data = t_half_model,
                      log_half_life ~  sigma_para_scaled +
                        NH_IR_scaled)

# Create a model that relates the sigma values to the log of the half life
lm_backbones_sigma <- lm(data = t_half_model,
                         log_half_life ~  sigma_para_scaled +
                           sigma_i_ar_scaled +
                           sigma_i_n_scaled)

# Create the same model but for logk instead of logthalf
lm_logk_backbones_sigma <- lm(data = t_half_model,
                              logk ~  sigma_para_scaled +
                                sigma_i_ar_scaled +
                                sigma_i_n_scaled)

#=========
# Plotting
#=========
# Plot showing concentrations over time for the various backbones
backbone_t_plot <- ggplot(data = half_life_data,
                          aes(x = time_s,
                              y = Pd_II_mol_norm)) +
  geom_point(aes(color = Backbone),
             size = 3) +
  scale_color_brewer(palette = "Paired") +
  labs(x = "Time (s)",
       y = "Relative Pd(II) Concentration") +
  themething

# Histogram of half lives and log half lives
t_half_hist <- ggplot(data = backbone_half_lives_PtBu3,
                      aes(x = half_life_s)) +
  geom_histogram(binwidth = 1e6) +
  geom_density(fill = "orange",
               alpha = 0.2) +
  labs(x = "Half Life at 25 \u00B0C (s)",
       y = "Count") +
  themething

log_t_half_hist <- ggplot(data = backbone_half_lives_PtBu3,
                          aes(x = log_half_life)) +
  geom_histogram(binwidth = 1) +
  geom_density(fill = "orange",
               alpha = 0.2) +
  labs(x = expression(paste(log[10]("Half-life"), " at 25 \u00B0C")),
       y = "Count") +
  themething

hist_comb <- t_half_hist +
  log_t_half_hist +
  plot_layout(axis_titles = "collect")

# Model summary plot for IR vs sigma model
model_comparison_plot <- plot_summs(lm_backbones_IR,
                                    lm_backbones_sigma,
                                    lm_logk_backbones_sigma,
                                    omit.coefs = NULL,
                                    plot.distributions = TRUE,
                                    model.names = c("N-H IR Stretch Model",
                                                    "Sigma I model",
                                                    "Sigma I model logk"),
                                    coefs = c("Intercept" = "(Intercept)",
                                              "Sigma para" = "sigma_para_scaled",
                                              "N-H IR Stretch" = "NH_IR_scaled",
                                              "Sigma I Ar" = "sigma_i_ar_scaled",
                                              "Sigma I N" = "sigma_i_n_scaled"))

model_comparison_table <- modelsummary(list(
  "IR Model" = lm_backbones_IR,
  "Sigma i model" = lm_backbones_sigma,
  "logk Model" = lm_logk_backbones_sigma),
fmt = 3,
estimate = "{estimate} [{conf.low}, {conf.high}]",
statistic = NULL,
gof_omit = "AIC|BIC|Log.Lik.",
title = "Summary Statistics for the Half Life Models",
output = "latex")

# Function for plotting the correlations from which I calculate the sigma_i
# values for the various substituents - needed for SI
# The arguments it takes are df_data - data frame your data is in
# x_var - whatever you want plotted on the x axis; y_var - cf. x_var
# pt_size - size of points, pt_col - color of points, x_lab - string 
# that is x axis label; y_lab - cf. x_lab
simple_plot_mb_style <- function(df_data,
                                 x_var,
                                 y_var,
                                 pt_size,
                                 pt_col,
                                 x_lab,
                                 y_lab)
{
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

# Plot for correlation between 19F NMR signals from Taft paper
taft_19F_plot <- simple_plot_mb_style(taft_ref_data,
                                      delta_F,
                                      sigma_i,
                                      3,
                                      "black",
                                      expression(
                                        paste({}^19,
                                              "F{",
                                              {}^1,
                                              "H} NMR signal relative to PhF")),
                                      expression(sigma[I]))

# Plot for correlation between pKa and sigma i
taft_pKa_plot <- simple_plot_mb_style(taft_pKa_data,
                                      pKa_H2O,
                                      sigma_i,
                                      3,
                                      "black",
                                      expression(paste(pK[a],
                                                       "(",
                                                       H[2],
                                                       "O)")),
                                      expression(sigma[I]))

# Plot for correlation between sigma i of sulfonyls and alkyls
sulfonyl_alkyl_sigma_plot <- simple_plot_mb_style(sigma_sulf,
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

# Plot backbone vs values
backbone_IR_plot <- simple_plot_mb_style(t_half_model,
                                         Backbone,
                                         N_H_stretch_IR,
                                         3,
                                         "black",
                                         "Backbone",
                                         expression(
                                           paste(
                                             "N-H IR Stretching Frequency (", 
                                             {cm}^{-1},
                                             ")"))) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1)) +
  scale_x_discrete(labels = parse(text = 
                                    t_half_model$Special_label_exp))

# Create a dataframe which contains predicted values and real values for 
# half life
t_half_pred_df <- t_half_model %>%
  mutate(log_t_half_pred = predict(lm_backbones_sigma,
                                   t_half_model))

# Create the plot for predicted vs measured half-life
t_half_pred_plot <- ggplot(data = t_half_pred_df,
                           aes(x = log_half_life,
                               y = log_t_half_pred)) +
  geom_abline(slope = 1,
              intercept = 0,
              color = "#4E525A",
              linewidth = 1) +
  geom_point(size = 5,
             color = "#AF0000") +
  geom_label_repel(data = filter(t_half_pred_df,
                                 (log_t_half_pred/log_half_life) >= 1
                                 | Backbone == "F4-NTf-OMe"),
                   aes(label = Special_label_exp),
                   nudge_x = -1,
                   nudge_y = 0.4,
                   parse = T,
                   size = 5) +
  geom_label_repel(data = filter(t_half_pred_df,
                                 (log_t_half_pred/log_half_life) < 1
                                 & Backbone != "F4-NTf-OMe"),
                   aes(label = Special_label_exp),
                   nudge_x = 1,
                   nudge_y = -0.4,
                   parse = T,
                   size = 5) +
  labs(x = expression(paste(log[10](t[1/2])," at 25 \u00B0C")),
       y = expression(paste("Predicted ", log[10](t[1/2])," at 25 \u00B0C"))) +
  xlim(0,8) +
  ylim(0,8) +
  themething +
  theme(axis.title = element_text(size = rel(1.75)),
        axis.text = element_text(size = rel(1.75)))