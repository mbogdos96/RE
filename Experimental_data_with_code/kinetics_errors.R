library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(ggbeeswarm)

source("SW_plots.R")
source("common_functions.R")
source("backbone_half_lives.R")
source("eyring.R")

#==================
#==================
# Backbone kinetics
#==================
#==================

#===================
# Data manipulations
#===================
# Import the data and then convert values to SI units
kinetics_revisions_df <- read.csv(
  file.path(
    "./half_life_data/revisions_kinetics_F4-NMs-Ph.csv")) %>%
  mutate(Temperature_K = Temperature_C + 273.15) %>%
  mutate(Internal_standard_mol = 
           (Internal_standard_mass_mg*0.001)/
           internal_standard_molecular_weight) %>%
  mutate(Pd_II_mol= Pd_II_integral*Internal_standard_mol) %>%
  group_by(Sample) %>%
  mutate(Pd_II_mol_norm = Pd_II_mol/max(Pd_II_mol)) %>%
  ungroup() 

# Import the data from the original run
all_F4NMsPh_df <- half_life_data %>%
  filter(Backbone == "F4-NMs-Ph") %>%
  mutate(Sample = "OG",
         time_s = time_s - 540) %>%
  dplyr::select(names(kinetics_revisions_df)) %>%
  rbind(kinetics_revisions_df)

# Get rate constants
k_r_rc_df <- calculate_rate_constant(
  filter(all_F4NMsPh_df,
         Sample == "S1"),
  0.1,
  1e-7) %>%
  add_row(calculate_rate_constant(filter(all_F4NMsPh_df,
                                         Sample == "S2"),
                                  0.1,
                                  1e-7)) %>%
  add_row(calculate_rate_constant(filter(all_F4NMsPh_df,
                                 Sample == "OG"),
                          0.1,
                          1e-7)) %>%
  mutate(Sample = c("S1","S2","OG"))

# Plot to show different datasets
kin_rep_pred_df <- data.frame(t = seq(0, range(
  filter(all_F4NMsPh_df, 
         Sample == "OG")$time_s)[[2]], 
  by = 10), 
  Sample = "OG") %>%
  rbind(data.frame(t = seq(0, 
                           range(
    filter(all_F4NMsPh_df, 
           Sample == "S1")$time_s)[[2]], 
    by = 10), 
    Sample = "S1"),
    data.frame(t = seq(0, 
                       range(
      filter(all_F4NMsPh_df, 
             Sample == "S2")$time_s)[[2]], 
      by = 10), 
      Sample = "S2")) %>%
  left_join(k_r_rc_df, 
            by = "Sample") %>%
  group_by(Sample) %>%
  mutate(Pd_pred = A0*exp(-rate_constant*t), 
         Pd_pred_norm = if_else(Sample == "S1",
 Pd_pred/max(
   data = filter(all_F4NMsPh_df,
                 Sample == "S1")$Pd_II_mol),
 if_else(Sample == "S2",
         Pd_pred/max(
           data = filter(all_F4NMsPh_df,
                         Sample == "S2")$Pd_II_mol),
         Pd_pred/max(
           data = filter(all_F4NMsPh_df,
                         Sample == "OG")$Pd_II_mol)))) %>%
  ungroup()

#=========
# Plotting
#=========
# Plot raw data and fit
kin_replicates_plot <- ggplot() +
  geom_line(data = kin_rep_pred_df,
            aes(x = t,
                y = Pd_pred_norm,
                color = Sample),
            linewidth = 2) +
  geom_point(data = all_F4NMsPh_df,
             aes(x = time_s,
                 y = Pd_II_mol_norm,
                 color = Sample),
             alpha = 0.5,
             size = 3) +
  scale_color_manual(values = c("#636C9D",
                                "#B71B1B",
                                "gray"),
                     labels = c("1",
                                "2",
                                "3"),
                     name = "Replicate") +
  labs(x = "Time (s)",
       y = "Relative [Pd(II)]") +
  themething

# Get mean and sd of k
k_err_df <- k_r_rc_df %>%
  mutate(mean_k = mean(rate_constant),
         upr_k = mean(rate_constant) 
         + sd(rate_constant),
         lwr_k = mean(rate_constant) -
           sd(rate_constant),
         upr_k_sem = mean(rate_constant) 
         + sd(rate_constant)/sqrt(3),
         lwr_k_sem = mean(rate_constant) 
         - sd(rate_constant)/sqrt(3))

# Raw observations and fit plot
mean_repl_plot <- ggplot() +
  geom_errorbar(data = k_err_df,
                aes(ymax = upr_k,
                    ymin = lwr_k,
                    x = Backbone),
                linewidth = 0.75,
                color = "black",
                width = 0.2) +
  geom_errorbar(data = k_err_df,
                aes(ymax = upr_k_sem,
                    ymin = lwr_k_sem,
                    x = Backbone),
                linewidth = 5,
                color = "blue",
                width = 0,
                alpha = 0.3) +
  geom_point(data = k_r_rc_df,
             aes(y = rate_constant,
                 x = Backbone,
                 color = Sample),
             size = 5) +
  scale_color_manual(values = c("#636C9D",
                                "#B71B1B",
                                "gray"),
                     labels = c("1",
                                "2",
                                "3"),
                     name = "Replicate") +
  ylim(0.8e-4,1.1e-4) +
  labs(y = expression(paste("k (",
                            s^{"-1"},
                            ")"))) +
  themething +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")

# Compare errors of logs for different backbones
logk_error_df <- k_err_df %>%
  mutate(logk_rt = log10(rate_constant),
         mean_logk = log10(mean_k),
         logk_upr = mean_logk + sd(logk_rt),
         logk_lwr = mean_logk - sd(logk_rt),
         logk_upr_sem = mean_logk + sd(logk_rt)/sqrt(3),
         logk_lwr_sem = mean_logk - sd(logk_rt)/sqrt(3),
         repeated = "Y")

logk_err_combined_df <- bind_rows(
  logk_error_df,
  backbone_half_lives_PtBu3 %>%
    mutate(logk_rt = log10(rate_constant))
)

# Plot of logs and error bars
logk_comp_plot <- ggplot() +
  geom_errorbar(data = logk_err_combined_df,
                aes(ymax = logk_upr,
                    ymin = logk_lwr,
                    x = repeated),
                linewidth = 0.75,
                color = "black",
                width = 0.5) +
  geom_errorbar(data = logk_err_combined_df,
                aes(ymax = logk_upr_sem,
                    ymin = logk_lwr_sem,
                    x = repeated),
                linewidth = 0.5,
                color = "blue",
                width = 1,
                alpha = 0.3) +
  geom_beeswarm(data = logk_err_combined_df,
             aes(y = logk_rt,
                 x = repeated,
                 color = Sample),
             size = 2) +
  scale_color_manual(values = c("#636C9D",
                                "#B71B1B",
                                "gray")) +
  themething +
  labs(y = expression(paste(log[10],
                            "(",
                            {}^{rt},
                            "k)"))) +
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        axis.ticks = element_blank(),
        axis.title.x = element_blank())

# Combine plots together
k_error_layout <- "AA
BC"

repl_kin_master_plot <- kin_replicates_plot +
  mean_repl_plot + 
  logk_comp_plot +
  plot_layout(nrow = 2,
              ncol = 1,
              design = k_error_layout)

#========================
#========================
# In situ kinetics PAd2Bu
#========================
#========================

#==============
# Data handling
#==============
# Import second set of data
repl_PAd2Bu <- read.csv(file.path(
  "./eyring/PAd2Bu_replicates.csv"
))

# Create list of 3 dfs one for each replicate
PAd2Bu_repl <-  repl_PAd2Bu %>%
  rename(time_s = "time...s",
         Pd_molar = "c.Pd.Complex....M") %>%
  mutate(Backbone = if_else(Replicate == 1,
                            "F4-NClMs-OMe1",
                            if_else(Replicate == 2,
                                    "F4-NClMs-OMe2",
                                    "F4-NClMs-OMe3"))) %>%
  group_by(Backbone) %>%
  mutate(rel_Pd_II = Pd_molar/max(Pd_molar))

# Splir into separate dfs
PAd2Bu_repl_list <- PAd2Bu_repl %>%
  group_split(Replicate)

# Apply fitting function used for other data
PAd2Bu_eyring <- fit_k_eyring(PAd2Bu_repl_list)

# Get mean and bounds for k
PAd2Bu_repl_df_k <- PAd2Bu_eyring$Rate_constants_eyring %>%
  mutate(avg_k = mean(rate_constant),
         upr_k = avg_k + sd(rate_constant),
         lwr_k = avg_k - sd(rate_constant),
         upr_k_sem = avg_k + sd(rate_constant)/sqrt(3),
         lwr_k_sem = avg_k - sd(rate_constant)/sqrt(3),
         logk_rt = log10(rate_constant),
         avg_logk = log10(avg_k),
         upr_logk = avg_logk + sd(logk_rt),
         lwr_logk = avg_logk - sd(logk_rt),
         logk_upr_sem = avg_logk + sd(logk_rt)/sqrt(3),
         logk_lwr_sem = avg_logk - sd(logk_rt)/sqrt(3),
         repeats = "Y",
         sample = c(1,2,3))

#=========
# Plotting
#=========
# Create df for plotting
eyring_comb <- bind_rows(
  PAd2Bu_repl_df_k %>%
    mutate(experimental = "Y"),
  eyring_summary %>%
    mutate(logk_rt = log10(rate_constant),
           repeats = "N",
           experimental = "Y"),
  k_E_L_scaled %>%
    mutate(experimental = "N",
           repeats = "N")) %>%
  mutate(experimental = as.factor(experimental),
         experimental = relevel(experimental,
                                ref = "Y"))

# Plot of repeats on their own
mean_repl_plot_eyring <- ggplot(data = eyring_comb %>%
                           filter(repeats == "Y")) +
  geom_errorbar(aes(ymax = upr_k,
                    ymin = lwr_k,
                    x = repeats),
                linewidth = 0.75,
                color = "black",
                width = 0.2) +
  geom_errorbar(aes(ymax = upr_k_sem,
                    ymin = lwr_k_sem,
                    x = repeats),
                linewidth = 5,
                color = "#636C9D",
                width = 0,
                alpha = 0.2) +
  geom_point(aes(y = rate_constant,
                 x = repeats,
                 colour = as.factor(sample)),
             size = 5) +
  labs(y = expression(paste("k (",
                            s^{"-1"},
                            ")"))) +
  scale_color_manual(values = c("#636C9D",
                                "#B71B1B",
                                "gray"),
                     labels = c("1",
                                "2",
                                "3"),
                     name = "Replicate") +
  themething +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())

# Plot to compare logks
logk_comp_plot_eyring <- ggplot(data = eyring_comb) +
  geom_errorbar(aes(ymax = upr_logk,
                    ymin = lwr_logk,
                    x = repeats),
                linewidth = 0.75,
                color = "black",
                width = 0.5) +
  geom_errorbar(aes(ymax = logk_upr_sem,
                    ymin = logk_lwr_sem,
                    x = repeats),
                linewidth = 0.5,
                color = "blue",
                width = 1,
                alpha = 0.3) +
  geom_beeswarm(aes(y = logk_rt,
                    x = repeats,
                    color = as.factor(sample),
                    shape = experimental),
                size = 2) +
  scale_color_manual(values = c("#636C9D",
                                "#B71B1B",
                                "gray")) +
  themething +
  labs(y = expression(paste(log[10],
                            "(k)"))) +
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        axis.ticks = element_blank(),
        axis.title.x = element_blank())

# Give relative concentrations for predictions
PAd2Bu_pred <- PAd2Bu_eyring$Predicted_values_eyring %>%
  group_by(Backbone) %>%
  mutate(Pd_rel = if_else(Backbone == "F4-NClMs-OMe1",
                          Pd_predicted/max(filter(PAd2Bu_repl,
                                               Backbone == 
                                                 "F4-NClMs-OMe1"
                                             )$Pd_molar),
                          if_else(Backbone == 
                                    "F4-NClMs-OMe2",
                                  Pd_predicted/max(
                                      filter(PAd2Bu_repl,
                                        Backbone == 
                            "F4-NClMs-OMe2")$Pd_molar),
                                  Pd_predicted/max(
                                      filter(PAd2Bu_repl,
                                        Backbone == 
                             "F4-NClMs-OMe3")$Pd_molar))))

# Plot of traces
kin_replicates_plot_eyring <- ggplot() +
  geom_line(data = PAd2Bu_pred,
            aes(x = time_s,
                y = Pd_rel,
                color = Backbone),
            linewidth = 2) +
  geom_point(data = PAd2Bu_repl,
             aes(x = time_s,
                 y = rel_Pd_II,
                 color = Backbone),
             alpha = 0.5,
             size = 3) +
  scale_color_manual(values = c("#636C9D",
                                "#B71B1B",
                                "gray"),
                     labels = c("1",
                                "2",
                                "3"),
                     name = "Replicate") +
  labs(x = "Time (s)",
       y = "Relative [Pd(II)]") +
  themething

# combine all plots
repl_kin_comb_plot_eyring <- kin_replicates_plot_eyring +
  mean_repl_plot_eyring + 
  logk_comp_plot_eyring +
  plot_layout(nrow = 2,
              ncol = 1,
              design = k_error_layout)