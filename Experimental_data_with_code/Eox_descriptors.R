library(dplyr)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(stringr)

source("common_functions.R")
source("SW_plots.R")

#====================================
#====================================
#====================================
# DFT descriptors vs Eox experimental
#====================================
#====================================
#====================================

#==============
#==============
# Data handling
#==============
#==============

#===================
# Free L descriptors
#===================

#============
# Descriptors
#============
# Keep only the oxidation potentials and the 
# ligand for simplicity
echem_data <- first_oxidation_potential_all %>%
  dplyr::select(Ligand,
                Potential)

# Import the CEP data from machine learning
CEP_ligands_ML <- read.csv(file.path(
  "./ligand_descriptors/ML_TEP.csv"))

# Import DFT calculated frequencies
TEP_DFT <- read.csv(file.path(
  "./ligand_descriptors/DFT_TEP.csv"))

# Create a dataframe which contains the DFT and ML TEPs
CO_stretch_method_df <- merge(CEP_ligands_ML,
                              TEP_DFT) %>%
  filter(L_type == "monodentate") %>%
  dplyr::select(Ligand,
                TolmanPrediction,
                Ni_CO_stretch) %>%
  mutate(Ni_CO_stretch = as.numeric(Ni_CO_stretch),
         rel_ML = 
           TolmanPrediction - min(TolmanPrediction),
         rel_DFT = 
           Ni_CO_stretch - min(Ni_CO_stretch),
         Special_label_exp = sapply(
           Ligand, 
           replace_ligand_expression_og))

# Create dataframe which contains both CO stretch and 
# oxidation potentials
Ni_Rh_TEP_E <- merge(TEP_DFT,
                     first_oxidation_potential_all) %>%
  mutate(Ni_CO_stretch = as.numeric(Ni_CO_stretch),
         Rh_CO_stretch = as.numeric(Rh_CO_stretch),
         rel_Ni_CO_stretch = Ni_CO_stretch - 
           min(Ni_CO_stretch,
               na.rm = T),
         rel_Rh_CO_stretch = Rh_CO_stretch 
         - min(Rh_CO_stretch,
               na.rm = T),
         Special_label_exp = sapply(
           Ligand, 
           replace_ligand_expression_og),
         rel_stretch = ifelse(
           !is.na(rel_Ni_CO_stretch),
           rel_Ni_CO_stretch,
           rel_Rh_CO_stretch))

# Check correlation CO stretched and palladacycle 
# oxidation potential for all, or each category
# Correlation for all
corr_TEP_E_all <- summary(lm(
  Potential ~ rel_stretch,
                             data = Ni_Rh_TEP_E))

# Correletion monodentate ligands
corr_TEP_E_mono <- summary(lm(
  Potential ~ Ni_CO_stretch,
                              data = Ni_Rh_TEP_E))

# Correlation bidentate ligands
corr_TEP_E_bi <- summary(lm(
  Potential ~ rel_Rh_CO_stretch,
                            data = Ni_Rh_TEP_E))

# Get free L oxidation potential together with
# palladacycle oxidation potential
pdcy_free_L_df <- first_ox_free_L %>%
  rename(Potential_free_L = "Potential",
         Current_free_L = "Current") %>%
  mutate(Ligand = str_replace(Ligand,
                              "_free$",
                              ""),
         Special_label_exp = sapply(
           Ligand, 
           replace_ligand_expression_og)) %>%
  left_join(first_oxidation_potential_all,
            by = "Ligand")

# Add data for DFT-calculated CO stretches
CO_EOx_free_L_df <- pdcy_free_L_df %>%
  left_join(Ni_Rh_TEP_E %>%
              dplyr::select(c(Ligand,
                              rel_Ni_CO_stretch,
                              rel_Rh_CO_stretch)),
            by = "Ligand")

# Import free L HOMO energies
ligand_HOMO_df <- read.csv(
  file.path(
    "./ligand_descriptors/ligand_HOMOs.csv")) %>%
  rename(HOMO_free_L_hartrees = "HOMO_energy_hartrees")

# Make one df for free L descriptors together
all_free_L_descr <- CO_EOx_free_L_df %>%
  left_join(ligand_HOMO_df %>%
              dplyr::select(-HOMO_energy_kcal_mol),
            by = "Ligand")

#========================
# Palldacycle descriptors
#========================
# Pd Mulliken charges and HOMO energies
Pd_charges <- read.csv(file.path(
  "./complex_descriptors/HOMO_and_charges.csv")) %>%
  merge(first_oxidation_potential_all %>%
          dplyr::select(Ligand,
                        Potential)) %>%
  mutate(HOMO_energy = as.numeric(HOMO_energy),
         Special_label_exp = sapply(
           Ligand, 
           replace_ligand_expression_og))

# Import VIP and merge with Eox data
VIP_desc_df <- read.csv(file.path(
  "./complex_descriptors/VIP_palldacycles.csv")) %>%
  dplyr::select(c(Ligand,
                  SP.VIP...eV,
                  AIP.DFT.CPCM....eV)) %>%
  left_join(Pd_charges,
            by = "Ligand")

#=========
#=========
# Plotting
#=========
#=========

#===================
# Free L descriptors
#===================
# Figure 4 plot for DFT TEP Ni or Rh vs Eox
E_v_TEP_Ni_Rh <- ggplot() +
  geom_label_repel(data = filter(Ni_Rh_TEP_E,
                  !is.na(rel_Ni_CO_stretch)),
                   aes(x = Potential,
                       y = rel_Ni_CO_stretch,
                       label = Special_label_exp),
                   nudge_x = 0.01,
                   nudge_y = 3,
                   parse = TRUE,
                   size = 4) +
  geom_label_repel(data = filter(Ni_Rh_TEP_E,
                                 !is.na(Rh_CO_stretch)),
                   aes(x = Potential,
                       y = rel_Rh_CO_stretch,
                       label = Special_label_exp),
                   nudge_x = 0.05,
                   nudge_y = -2,
                   parse = TRUE,
                   size = 4) +
  geom_point(data = filter(Ni_Rh_TEP_E,
                           !is.na(rel_Rh_CO_stretch)),
             aes(x = Potential,
                 y = rel_Rh_CO_stretch),
             fill = "#B71B1B",
             size = 4,
             pch = 21) +
  geom_point(data = filter(Ni_Rh_TEP_E,
                           !is.na(rel_Ni_CO_stretch)),
             aes(x = Potential,
                 y = rel_Ni_CO_stretch),
             fill = "#636C9D",
             size = 4,
             pch = 21) +
  labs(x = expression(paste("Palladacycle ", 
                            E^{ox} ," vs. Fc/", 
                            {"Fc"}^{"+"}, " (V)")),
       y = expression(paste("Rh(I) ", 
                            {v}[CO], " (rel., ",
                            "DFT, ", {cm}^{-1}, ")"))) +
  scale_y_continuous(sec.axis = sec_axis(~., 
                  name = expression(paste(
                    "Ni(0) ", {v}[CO], " (rel.",
                    ", DFT, ", {cm}^{-1}, ")")))) +
  coord_cartesian(xlim = c(0.7,1.3)) +
  themething +
  theme(axis.title.y.right = element_text(
    color = "#636C9D"),
        axis.title.y.left = element_text(
          color = "#B71B1B"))

# Plot comparing DFT vs ML TEP
TEP_DFT_ML_plot <- ggplot(data = CO_stretch_method_df) +
  geom_label_repel(aes(x = Ni_CO_stretch,
                       y = TolmanPrediction,
                       label = Special_label_exp),
                   parse = T,
                   nudge_x = 0.5,
                   nudge_y = 1,
                   size = 5) +
  labs(x = expression(paste(
    "DFT CO Stretching Frequency (", 
                            {cm}^{-1}, 
                            ")")),
       y = expression(paste(
         "Machine Learning CO Stretching Frequency (", 
                            {cm}^{-1}, 
                            ")"))) +
  geom_point(aes(x = Ni_CO_stretch,
                 y = TolmanPrediction),
             size = 5) +
  themething

# Plot comparing DFT vs ML TEP relative values
TEP_DFT_ML_plot_rel <- ggplot(data = 
                                CO_stretch_method_df) +
  geom_label_repel(aes(x = rel_DFT,
                       y = rel_ML,
                       label = Special_label_exp),
                   nudge_x = 0.05,
                   nudge_y = 0.1,
                   size = 5) +
  labs(x = expression(paste(
    "Relative DFT CO Stretching Frequency (", 
                            {cm}^{-1}, ")")),
       y = expression(paste(
"Relative Machine Learning CO Stretching Frequency (", 
         {cm}^{-1}, ")"))) +
  geom_point(aes(x = rel_DFT,
                 y = rel_ML),
             size = 5) +
  themething

# Plot of SW peak for free L and pdcy w/that L
# Function for plots
mk_desc_corr_plot <- function(inp_df,
                               aes_x,
                               aes_y,
                               nudge_x_1,
                               nudge_y_1,
                               nudge_x_2,
                               nudge_y_2,
                               lab_cond_1,
                               lab_cond_2,
                               add_correlation,
                               add_xeqy){
  
  # Initiate empty plot
  out_plot <- ggplot(data = inp_df,
                     aes(x = {{aes_x}},
                         y = {{aes_y}}))
  
  # Optionally add a correlation
  if(isTRUE(add_correlation)){
    out_plot <- out_plot +
      geom_smooth(method = lm,
                  se = T,
                  alpha = 0.3,
                  color = "black")
  }
  
  # Otherwise do not modify plot
  else{
    out_plot <- out_plot
  }
  
  # Optionally add x = y line
  if(isTRUE(add_xeqy)){
    # Add other line
    out_plot <- out_plot +
      geom_abline(slope = 1,
                  intercept = 0,
                  color = "#4E525A",
                  linewidth = 1)
  }
  
  # Otherwise do not modify the plot
  else{
    out_plot <- out_plot
  }
  
  # Modify labels differentially
  # If you choose false the label will be
  # applied the same to all datapoints
  if(isFALSE(lab_cond_1)){
    out_plot <- out_plot +
      geom_label_repel(aes(label = Special_label_exp),
                       nudge_x = nudge_x_1,
                       nudge_y = nudge_y_1,
                       parse = T,
                       size = 4)
  }
  
  # If you input different conditions then the labels
  # will be nudged differently for the two groups
  else{
    out_plot <- out_plot +
      geom_label_repel(data = {{lab_cond_1}},
                       aes(label = Special_label_exp),
                       nudge_x = nudge_x_1,
                       nudge_y = nudge_y_1,
                       parse = T,
                       size = 4) +
      geom_label_repel(data = {{lab_cond_2}},
                       aes(x = {{aes_x}},
                           y = {{aes_y}},
                           label = Special_label_exp),
                       nudge_x = nudge_x_2,
                       nudge_y = nudge_y_2,
                       parse = T, 
                       size = 4)
  }
  
  # Add other features of the plot
  out_plot <- out_plot +
    geom_point(aes(fill = L_type),
               size = 4,
               pch = 21) +
    scale_fill_manual(values = c("#B71B1B",
                                  "#636C9D")) +
    themething
  
  # Output the plot
  return(out_plot)
}

# Create the plot
free_L_complex_E_plot <- mk_desc_corr_plot(pdcy_free_L_df,
                                           Potential_free_L,
                                           Potential,
                                           0.1,
                                           0,
                                           0.1,
                                           0,
                                           F,F,F,F) +
  labs(x = expression(paste(
    E^{ox},
    " Free L (Fc/",
    {Fc}^{"+"},
    ", V)"
  )),
  y = expression(paste(
    E^{ox},
    " Palladacycle (Fc/",
    {Fc}^{"+"},
    ", V)"
  )),
  color = "Ligand Denticity") +
  themething +
  theme(legend.position = "none")

# Plot Eox vs CO stretches
plot_Eox_CO_free_L <- ggplot(data = CO_EOx_free_L_df,
                             aes(x = Potential_free_L)) +
  geom_label_repel(aes(y = rel_Rh_CO_stretch,
                       label = Special_label_exp),
                   nudge_y = 10,
                   parse = T,
                   size = 4) +
  geom_label_repel(aes(y = rel_Ni_CO_stretch,
                       label = Special_label_exp),
                   nudge_y = -10,
                   parse = T,
                   size = 4) +
  geom_point(aes(y = rel_Rh_CO_stretch,
                 fill = L_type),
             size = 4,
             pch = 21) +
  geom_point(aes(y = rel_Ni_CO_stretch,
                 fill = L_type),
             size = 4,
             pch = 21) +
  scale_fill_manual(values = c("#B71B1B",
                                "#636C9D")) +
  labs(x = expression(paste(
    E^{ox},
    " Free L (Fc/",
    {Fc}^{"+"},
    ", V)"
  )),
  x = expression(paste(
    E^{ox},
    " Palladacycle (Fc/",
    {Fc}^{"+"},
    ", V)"
  )),
  y = expression(paste("Rh(I) ", 
                       {v}[CO], " (rel., ",
                       "DFT, ", {cm}^{-1}, ")")),
  color = "Ligand Denticity") +
  scale_y_continuous(sec.axis = sec_axis(~., 
  name = expression(
    paste(
      "Ni(0) ", {v}[CO], " (rel.",
      ", DFT, ", {cm}^{-1}, ")")))) +
  themething +
  theme(legend.position = "none")

# Combined plot for Figure 5
E_ox_no_corr_plot <- plot_Eox_CO_free_L +
  free_L_complex_E_plot +
  plot_layout(nrow = 1,
              ncol = 2,
              axes = "collect")

# Correlation between free L HOMO and free L Eox
HOMO_Eox_free_L_plot <- mk_desc_corr_plot(
  all_free_L_descr,
  HOMO_free_L_hartrees,
  Potential_free_L,
  -0.003,
  -0.1,
  0.01,
  0.1,
  . %>%
    filter(
      L_type == "bidentate"),
  . %>%
    filter(
      L_type == "monodentate"),
  F,F) +
  theme(legend.position = "none") +
  labs(y = expression(paste("Free L Potential (V, Fc/", 
                            {Fc^{"+"}} ,
                            ")")),
       x = "Free L HOMO Energy (Hartrees)")

# Correlation between v(CO) and free L HOMO
plot_HOMO_CO <- ggplot(data = all_free_L_descr,
                             aes(x = HOMO_free_L_hartrees)) +
  geom_label_repel(aes(y = rel_Rh_CO_stretch,
                       label = Special_label_exp),
                   nudge_y = -5,
                   nudge_x = -0.005,
                   parse = T) +
  geom_label_repel(aes(y = rel_Ni_CO_stretch,
                       label = Special_label_exp),
                   nudge_y = 10,
                   nudge_x = 0.005,
                   parse = T) +
  geom_point(aes(y = rel_Rh_CO_stretch,
                 fill = L_type),
             size = 3,
             pch = 21) +
  geom_point(aes(y = rel_Ni_CO_stretch,
                 fill = L_type),
             size = 3,
             pch = 21) +
  scale_fill_manual(values = c("#B71B1B",
                               "#636C9D")) +
  labs(x = "Free L HOMO Energy (Hartrees)",
  y = expression(paste("Rh(I) ", 
                       {v}[CO], " (rel., ",
                       "DFT, ", {cm}^{-1}, ")")),
  color = "Ligand Denticity") +
  scale_y_continuous(sec.axis = sec_axis(~., 
                                         name = expression(
                                           paste(
                                             "Ni(0) ", {v}[CO], " (rel.",
                                             ", DFT, ", {cm}^{-1}, ")")))) +
  themething +
  theme(legend.position = "none")

# R1 HOMO plot
r1_homo_plot <- plot_HOMO_CO +
  HOMO_Eox_free_L_plot +
  plot_layout(axes = "collect")

#========================
# Palldacycle descriptors
#========================
# Plot Mulliken charges vs experimental Eox
Mull_plot <- mk_desc_corr_plot(Pd_charges,
                               Potential,
                               Pd_Mulliken_charge,
                               0.01,
                               -0.1,
                               -0.01,
                               0.1,
                               . %>%
                                 filter(
                                   L_type == "bidentate"),
                               . %>%
                                 filter(
                                   L_type == "monodentate"),
                               F,F) +
  labs(x = expression(paste(
    E^{ox},
    " Palladacycle (Fc/",
    {Fc}^{"+"},
    ", V)")),
       y = "Mulliken Charge at Pd") +
  themething +
  theme(legend.position = "none")

# HOMO energy vs Eox
HOMO_plot <- mk_desc_corr_plot(Pd_charges,
                               Potential,
                               HOMO_energy,
                               -0.01,
                               -0.002,
                               0.01,
                               0.002,
                               . %>%
                                 filter(
                                   L_type == "bidentate"),
                               . %>%
                                 filter(
                                   L_type == "monodentate"),
                               F,F) +
  labs(x = expression(paste(
    E^{ox},
    " Palladacycle (Fc/",
    {Fc}^{"+"},
    ", V)")),
       y = "HOMO Energy (Hartrees)") +
  themething +
  theme(legend.position = "none")

# Plot VIP vs Eox
VIP_plot <- mk_desc_corr_plot(VIP_desc_df,
                              Potential,
                              SP.VIP...eV,
                              -0.01,
                              -0.03,
                              0,
                              0.03,
                              . %>%
                                filter(
                                  L_type == "bidentate"),
                              . %>%
                                filter(
                                  L_type == "monodentate"),
                              F,F) +
  labs(x = expression(paste(
    E^{ox},
    " Palladacycle (Fc/",
    {Fc}^{"+"},
    ", V)")),
       y = "Vertical Ionization Potential (eV)") +
  themething +
  theme(legend.position = "none")

# AIP plot
AIP_plot <- mk_desc_corr_plot(VIP_desc_df,
                              Potential,
                              AIP.DFT.CPCM....eV,
                              -0.01,
                              -0.002,
                              0.01,
                              0.002,
                              . %>%
                                filter(
                                  L_type == "bidentate"),
                              . %>%
                                filter(
                                  L_type == "monodentate"),
                              F,F) +
  labs(x = expression(paste(
    E^{ox},
    " Palladacycle (Fc/",
    {Fc}^{"+"},
    ", V)")),
       y = "Adiabatic Ionization Potential (eV)") +
  themething +
  theme(legend.position = "none")

# Figure 12
desc_E_corr_plot <- Mull_plot +
  HOMO_plot + 
  VIP_plot +
  AIP_plot +
  plot_layout() &
  theme(axis.title = element_text(size = rel(0.9)),
        axis.text = element_text(size = rel(0.9)))
