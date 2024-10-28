library(dplyr)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(stringr)

source("SW_analysis.R")
source("common_functions.R")

#=======================
# Functions for plotting
#=======================
# Function for plotting the square wave data, 
# current vs potential. Requires:
# df = dataframe with square wave data
plot_SW <- function(df){
  SWs <- ggplot(data = df, 
                aes(x = Potential, 
                    y = Current, 
                    color = Ligand)) +
    geom_point() +
    labs(x = expression(paste("Potential (V, Fc/", 
                              {Fc^{"+"}} ,
                              ")")),
         y = "Current (A)", 
         color = "Ligand") +
    themething
}

# Function for plotting the values of 
# extracted peaks in square wave data
SW_peak_plot <- function(df){
  
  first_ox_plot <- ggplot(data = df,
                          aes(x = factor(Ligand,
                        levels = unique(Ligand)),
                              y = Potential,
                              fill = L_type)) +
    geom_point(shape = 22, 
               size = 7) +
    labs(x = "Ligand", 
         y = expression(paste("First ",
                              E["ox"],
                              " (V, Fc/", 
                              {Fc^{"+"}} ,
                              ")"))) +
    themething +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1))
}

#==============
# Various plots
#==============
# Square wave voltammetry plots
# All ligands
SW_plot_all <- plot_SW(all_SW)

# Separate bidentate and monodentate ligands 
# into two plots
SW_L_type <- (plot_SW(all_SW %>%
                        filter(L_type == "bidentate")) +
                xlim(0.4,1.75) +
                ylim(-5e-7,1e-5) +
                theme(legend.position = "bottom") +
                guides(color = guide_legend(nrow = 2))) + 
  (plot_SW(all_SW %>%
            filter(L_type == "monodentate")) +
     xlim(0.5,1.75) +
     ylim(-5e-7,1e-5) +
  theme(legend.position = "bottom",
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
  guides(color = guide_legend(nrow = 2))) +
  plot_layout(axis_titles = "collect")

# Plots of square wave extracted peaks
first_oxidation_potential_all_mod <- 
  first_oxidation_potential_all %>%
  mutate(Special_label_exp = sapply(
    Ligand,
    replace_ligand_expression_og))

# First oxidation potential
E1_plot_all <- SW_peak_plot(
  first_oxidation_potential_all_mod) +
  scale_x_discrete(labels = parse(
    text = 
    first_oxidation_potential_all_mod$Special_label_exp)) +
  labs(y = expression(paste(E["ox"],
                            " (V, Fc/", 
                            {Fc^{"+"}} ,
                            ")")))

# All oxidation potentials of peaks
first_three_oxidation_potentials_all_mod <- 
  first_three_oxidation_potentials_all %>%
  mutate(Special_label_exp = sapply(
    Ligand, 
    replace_ligand_expression_og))

# Make the plot
E3_plot_all <- SW_peak_plot(
  first_three_oxidation_potentials_all_mod) +
  scale_x_discrete(labels = parse(
    text = 
    first_oxidation_potential_all_mod$Special_label_exp)) +
  labs(y = expression(paste(E["ox"],
                            " (V, Fc/", 
                            {Fc^{"+"}} ,
                            ")")))
# Combine the plots
comb_E1_E3 <- E3_plot_all +
  E1_plot_all +
  plot_layout(guides = "collect",
              axes = "collect") &
  theme(legend.position = "none")

# SW free L monodentate
SW_free_L_mono <- plot_SW(all_free_L_SW %>%
                    mutate(Ligand = str_replace(Ligand,
                                                "_free$",
                                                "")) %>%
                        filter(Ligand %in% c("PtBu3",
                                             "PtBu2Np",
                                             "PAd2Bu",
                                             "PMes3",
                                             "PtBu2Cy"))) +
                xlim(0,1.6) +
                ylim(-5e-7,1.5e-5) +
                theme(legend.position = "bottom") +
                guides(color = guide_legend(nrow = 2))

# SW free L bidentate
SW_free_L_bid <- plot_SW(all_free_L_SW %>%
                            mutate(Ligand = str_replace(Ligand,
                                                        "_free$",
                                                        "")) %>%
                            filter(!(Ligand %in% c("PtBu3",
                                                 "PtBu2Np",
                                                 "PAd2Bu",
                                                 "PMes3",
                                                 "PtBu2Cy")))) +
  xlim(0,1.6) +
  ylim(-5e-7,1.5e-5) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2))

# SW free L combined plots
SW_free_L_plot <- SW_free_L_mono +
  SW_free_L_bid +
  plot_layout(axes = "collect")
