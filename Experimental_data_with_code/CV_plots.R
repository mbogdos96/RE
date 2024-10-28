library(ggplot2)
library(dplyr)
library(ggrastr)
library(patchwork)

source("CV_analysis.R")
source("SW_plots.R")
source("common_functions.R")

#==================
#==================
# Scan rate studies
#==================
#==================

#=================================================
# Plots current vs voltage at different scan rates
#=================================================
# Function for plotting the CV data with varying scan rates - Warning VERY SLOW!
plot_CV_scanrate <- function(df){
  plotlist <- list()
  
  ligands <- unique(df$Ligand)
  
  for (ligand in ligands) {
    f_df <- df %>%
      filter(Ligand == ligand) %>%
      filter(Scan == 2)
    
    CVs_scanrate <- ggplot(data = f_df, 
                           aes(x = Potential, 
                               y = Current, 
                               color = Scan_rate)) +
      rasterize(geom_point(size = 0.25),
                dpi = 75) +
      labs(x = expression(paste("Voltage (V, vs. Fc/", {Fc}^{"+"} ,")")), 
           y = "Current (A)") +
      themething
    
    plotlist[[as.character(ligand)]] <- CVs_scanrate
  }
  return(plotlist)
}

# Generate the plot for the isolated peaks
CVs_scanrate_plot_isolated_peaks <- plot_CV_scanrate(
  CVs_isolated_peaks_scanrate_df)

# Generate the plot for all peaks
CVs_scanrate_plot_all_peaks <- plot_CV_scanrate(
  CVs_all_peaks_scanrate_df)

#================================================
# Plots voltage of peak current vs log(scan rate)
#================================================
# Function which creates the graph of Potential vs. log(scan rate)
# for any number of ligands in the analysis. 
# Arguments:
# df1 - contains the values for potential and log(scan rate)
# value - list of colours to be used; should be same number and number of ligands
CV_scanrate_fit_n_plot <- function(df1,value_x,value_y) {
  # Create an empty ggplot object
  plot <- ggplot()
  
  # Iterate over unique ligands
  for (ligand in unique(df1$Ligand)) {
    # Subset dataframes for the current ligand
    group_df1 <- filter(df1, 
                        Ligand == ligand)
    
    # Add points to the plot
    plot <- plot +
      geom_point(data = group_df1, 
                 aes(x = {{value_x}},
                     y = {{value_y}},
                     color = Ligand))
    
    # Add smoothed line to the plot
    plot <- plot +
      geom_smooth(data = group_df1, 
                  method = "lm", 
                  se = FALSE, 
                  aes(x = {{value_x}}, 
                      y = {{value_y}}, 
                      color = Ligand))
  }
  
  # Apply additional plot settings
  plot <- plot +
    scale_color_brewer(palette = "Dark2") +
    themething
  
  return(plot)
}

# Identify ligands as mono or bidentate
l_scan_rate <- ligands_scanrate %>%
  mutate(L_type = ifelse(
    Ligand %in% c("BIPHEP", 
                  "dcype", 
                  "dppe", 
                  "Fdppe",
                  "dppm"
    ),
    "bidentate",
    "monodentate"))

# Plot of current against squate root scanrate
CV_sqrtSR_I_plot <- (CV_scanrate_fit_n_plot(filter(l_scan_rate,
                                                   L_type == "bidentate"),
                                            Scan_rate_sqrt,
                                            Current) +
                       ylim(0,6e-5) +
                       guides(color = guide_legend(nrow = 2))) +
  (CV_scanrate_fit_n_plot(filter(l_scan_rate,
                                 L_type == "monodentate"),
                          Scan_rate_sqrt,
                          Current) +
     ylim(0,6e-5) +
     theme(axis.ticks.y = element_blank(),
           axis.text.y = element_blank()) +
     guides(color = guide_legend(nrow = 2))) +
  plot_layout(axis_titles = "collect",
              ncol = 2,
              nrow = 1) &
  theme(legend.position = "bottom") & 
  labs(x = expression(paste(sqrt("Scan Rate"), 
                            " (V/s", 
                            {")"}^{0.5})),
       y = "Peak Current (A)")

# Plot of current against scan rate
CV_SR_I_plot <- (CV_scanrate_fit_n_plot(filter(l_scan_rate,
                                               L_type == "bidentate"),
                                        Scan_rate,
                                        Current) +
                   ylim(0,6e-5) +
                   guides(color = guide_legend(nrow = 2))) +
  (CV_scanrate_fit_n_plot(filter(l_scan_rate,
                                 L_type == "monodentate"),
                          Scan_rate,
                          Current) +
     ylim(0,6e-5) +
     theme(axis.ticks.y = element_blank(),
           axis.text.y = element_blank()) +
     guides(color = guide_legend(nrow = 2))) +
  plot_layout(axis_titles = "collect") &
  theme(legend.position = "bottom") &
  labs(x = "Scan Rate (V/s)",
       y = "Peak Current (A)")

# Plot of current against log10 scan rate
CV_logSR_V_plot <- (CV_scanrate_fit_n_plot(filter(l_scan_rate,
                                                  L_type == "bidentate"),
                                           Scan_rate_log,
                                           Potential) +
                      labs(x = expression(paste(log[10]("Scan Rate"), 
                                                " (",
                                                log[10]("V/s"),
                                                ")")),
                           y = expression(paste("First ",
                                                E["ox"],
                                                " (V, Fc/", 
                                                {Fc^{"+"}} ,
                                                ")"))) +
                      ylim(0.65,1.5) +
                      theme(legend.position = "bottom") +
                      guides(color = guide_legend(nrow = 2))) +
  (CV_scanrate_fit_n_plot(filter(l_scan_rate,
                                 L_type == "monodentate"),
                          Scan_rate_log,
                          Potential) +
     labs(x = expression(paste(log[10]("Scan Rate"), 
                               " (",
                               log[10]("V/s"),
                               ")"))) +
     ylim(0.65,1.5) +
     theme(legend.position = "bottom",
           axis.ticks.y = element_blank(),
           axis.text.y = element_blank(),
           axis.title.y = element_blank()) +
     guides(color = guide_legend(nrow = 2))) +
  plot_layout(axis_titles = "collect")

#=======================
#=======================
# Scan direction studies
#=======================
#=======================
# Function for plotting the CV data with varying scan directions
plot_CV_scandir <- function(df){
  plotlist <- list()
  
  ligands <- unique(df$Ligand)
  
  for (ligand in ligands) {
    f_df <- df %>%
      filter(Ligand == ligand) %>%
      filter(Scan == 1)
    
    CVs_scandir <- ggplot(data = f_df, 
                          aes(x = Potential, 
                              y = Current, 
                              color = Scan_direction)) +
      rasterize(geom_point(size = 0.25),
                dpi = 75) +
      labs(x = expression(paste("Voltage (V, vs. Fc/", {Fc}^{"+"}," )")), 
           y = "Current (A)") +
      themething +
      theme(legend.position = "none")
    
    plotlist[[as.character(ligand)]] <- CVs_scandir
  }
  
  return(plotlist)
}

# Generate the plot for the scan directions
CVs_dir_plot <- plot_CV_scandir(CV_dir)

# Function for just showing the positive second scan
plot_CV_pos_2 <- function(df){
  plotlist <- list()
  
  ligands <- unique(df$Ligand)
  
  for (ligand in ligands) {
    f_df <- df %>%
      filter(Ligand == ligand,
             Scan == 2,
             Scan_direction == "pos")
    
    CVs_scandir <- ggplot(data = f_df, 
                          aes(x = Potential, 
                              y = Current)) +
      rasterize(geom_point(size = 0.25),
                dpi = 75) +
      labs(x = expression(paste("Voltage (V, vs. Fc/", 
                                {Fc}^{"+"}," )")), 
           y = "Current (A)") +
      themething +
      theme(legend.position = "none")
    
    plotlist[[as.character(ligand)]] <- CVs_scandir
  }
  
  return(plotlist)
}

# Generate the plots for the second scan of the CV
# done in the positive direction
CVs_pub <- plot_CV_pos_2(CV_dir)

# Combine plots of ligands to make them nice for the SI
# Nothing interesting from here on down

# dcype
# Regular CV
CV_dcype <- CVs_pub$dcype +
  annotate("segment",
           x = -0.3,
           y = 2e-06,
           xend = 0,
           yend = 2e-06,
           linewidth = 2,
           arrow = arrow(length = unit(0.15,
                                       "inches")))

# Scan direction
CV_dcype_dir <- CVs_dir_plot$dcype +
  annotate("segment",
           x = -0.3,
           y = 2e-06,
           xend = 0,
           yend = 2e-06,
           linewidth = 2,
           arrow = arrow(length = unit(0.15,
                                       "inches")),
           color = "#00BEC4") +
  annotate("segment",
           x = 0,
           y = 3e-06,
           xend = -0.3,
           yend = 3e-06,
           linewidth = 2,
           arrow = arrow(length = unit(0.15,
                                       "inches")),
           color = "#F8766D")

# Combine with patchwork
CV_studies_dcype <- CV_dcype +
  CV_dcype_dir +
  plot_layout(axis_titles = "collect")

# dppe
# Regular CV
CV_dppe <- CVs_pub$dppe +
  annotate("segment",
           x = -0.3,
           y = 0.5e-05,
           xend = 0,
           yend = 0.5e-05,
           linewidth = 2,
           arrow = arrow(length = unit(0.15,
                                       "inches")))

# Scan direction
CV_dppe_dir <- CVs_dir_plot$dppe +
  annotate("segment",
           x = -0.3,
           y = 0.5e-05,
           xend = 0,
           yend =  0.5e-05,
           linewidth = 2,
           arrow = arrow(length = unit(0.15,
                                       "inches")),
           color = "#00BEC4") +
  annotate("segment",
           x = 0,
           y = 0.75e-05,
           xend = -0.3,
           yend = 0.75e-05,
           linewidth = 2,
           arrow = arrow(length = unit(0.15,
                                       "inches")),
           color = "#F8766D")

# Combine with patchwork
CV_studies_dppe <- CV_dppe +
  CV_dppe_dir +
  plot_layout(axis_titles = "collect")

# dppe
# Regular CV
CV_Fdppe <- CVs_pub$Fdppe +
  annotate("segment",
           x = -0.3,
           y = 0.5e-05,
           xend = 0,
           yend = 0.5e-05,
           linewidth = 2,
           arrow = arrow(length = unit(0.15,
                                       "inches")))

# Scan direction
CV_Fdppe_dir <- CVs_dir_plot$Fdppe +
  annotate("segment",
           x = -0.3,
           y = 0.5e-05,
           xend = 0,
           yend =  0.5e-05,
           linewidth = 2,
           arrow = arrow(length = unit(0.15,
                                       "inches")),
           color = "#00BEC4") +
  annotate("segment",
           x = 0,
           y = 0.75e-05,
           xend = -0.3,
           yend = 0.75e-05,
           linewidth = 2,
           arrow = arrow(length = unit(0.15,
                                       "inches")),
           color = "#F8766D")

# Combine with patchwork
CV_studies_Fdppe <- CV_Fdppe +
  CV_Fdppe_dir +
  plot_layout(axis_titles = "collect")

# BIPHEP
# Regular CV
CV_BIPHEP <- CVs_pub$BIPHEP +
  annotate("segment",
           x = -0.3,
           y = 0.5e-05,
           xend = 0,
           yend = 0.5e-05,
           linewidth = 2,
           arrow = arrow(length = unit(0.15,
                                       "inches")))
# Scan direction
CV_BIPHEP_dir <- CVs_dir_plot$BIPHEP +
  annotate("segment",
           x = -0.3,
           y = 0.5e-05,
           xend = 0,
           yend =  0.5e-05,
           linewidth = 2,
           arrow = arrow(length = unit(0.15,
                                       "inches")),
           color = "#00BEC4") +
  annotate("segment",
           x = 0,
           y = 0.75e-05,
           xend = -0.3,
           yend = 0.75e-05,
           linewidth = 2,
           arrow = arrow(length = unit(0.15,
                                       "inches")),
           color = "#F8766D")

# Combine with patchwork
CV_studies_BIHPEP <- CV_BIPHEP +
  CV_BIPHEP_dir +
  plot_layout(axis_titles = "collect")

# dppm
# Regular CV
CV_dppm <- CVs_pub$dppm +
  annotate("segment",
           x = -0.3,
           y = 3e-06,
           xend = 0,
           yend = 3e-06,
           linewidth = 2,
           arrow = arrow(length = unit(0.15,
                                       "inches")))

# Scan direction
CV_dppm_dir <- CVs_dir_plot$dppm +
  annotate("segment",
           x = -0.3,
           y = 3e-06,
           xend = 0,
           yend = 3e-06,
           linewidth = 2,
           arrow = arrow(length = unit(0.15,
                                       "inches")),
           color = "#00BEC4") +
  annotate("segment",
           x = 0,
           y = 4e-06,
           xend = -0.3,
           yend = 4e-06,
           linewidth = 2,
           arrow = arrow(length = unit(0.15,
                                       "inches")),
           color = "#F8766D")

# Combine with patchwork
CV_studies_dppm <- CV_dppm +
  CV_dppm_dir +
  plot_layout(axis_titles = "collect")

# PAd2Bu
# Regular CV
CV_PAd2Bu <- CVs_pub$PAd2Bu +
  annotate("segment",
           x = -0.3,
           y = 3e-06,
           xend = 0,
           yend = 3e-06,
           linewidth = 2,
           arrow = arrow(length = unit(0.15,
                                       "inches")))

# Scan direction
CV_PAd2Bu_dir <- CVs_dir_plot$PAd2Bu +
  annotate("segment",
           x = -0.3,
           y = 3e-06,
           xend = 0,
           yend = 3e-06,
           linewidth = 2,
           arrow = arrow(length = unit(0.15,
                                       "inches")),
           color = "#00BEC4") +
  annotate("segment",
           x = 0,
           y = 4e-06,
           xend = -0.3,
           yend = 4e-06,
           linewidth = 2,
           arrow = arrow(length = unit(0.15,
                                       "inches")),
           color = "#F8766D")

# Combine with patchwork
CV_studies_PAd2Bu <- CV_PAd2Bu +
  CV_PAd2Bu_dir +
  plot_layout(axis_titles = "collect")

# PMes3
# Regular CV
CV_PMes3 <- CVs_pub$PMes3 +
  annotate("segment",
           x = -0.3,
           y = 4e-06,
           xend = 0,
           yend = 4e-06,
           linewidth = 2,
           arrow = arrow(length = unit(0.15,
                                       "inches")))

# Scan direction
CV_PMes3_dir <- CVs_dir_plot$PMes3 +
  annotate("segment",
           x = -0.3,
           y = 4e-06,
           xend = 0,
           yend = 4e-06,
           linewidth = 2,
           arrow = arrow(length = unit(0.15,
                                       "inches")),
           color = "#00BEC4") +
  annotate("segment",
           x = 0,
           y = 6.5e-06,
           xend = -0.3,
           yend = 6.5e-06,
           linewidth = 2,
           arrow = arrow(length = unit(0.15,
                                       "inches")),
           color = "#F8766D")

# Combine with patchwork
CV_studies_PMes3 <- CV_PMes3 +
  CV_PMes3_dir +
  plot_layout(axis_titles = "collect")


# PtBu2Cy
# Regular CV
CV_PtBu2Cy <- CVs_pub$PtBu2Cy +
  annotate("segment",
           x = -0.3,
           y = 2.5e-06,
           xend = 0,
           yend = 2.5e-06,
           linewidth = 2,
           arrow = arrow(length = unit(0.15,
                                       "inches")))

# Scan direction
CV_PtBu2Cy_dir <- CVs_dir_plot$PtBu2Cy +
  annotate("segment",
           x = -0.3,
           y = 2.5e-06,
           xend = 0,
           yend = 2.5e-06,
           linewidth = 2,
           arrow = arrow(length = unit(0.15,
                                       "inches")),
           color = "#00BEC4") +
  annotate("segment",
           x = 0,
           y = 4e-06,
           xend = -0.3,
           yend = 4e-06,
           linewidth = 2,
           arrow = arrow(length = unit(0.15,
                                       "inches")),
           color = "#F8766D")

# Combine with patchwork
CV_studies_PtBu2Cy <- CV_PtBu2Cy +
  CV_PtBu2Cy_dir +
  plot_layout(axis_titles = "collect")


# PtBu2Np
# Regular CV
CV_PtBu2Np <- CVs_pub$PtBu2Np +
  annotate("segment",
           x = -0.3,
           y = 2.5e-06,
           xend = 0,
           yend = 2.5e-06,
           linewidth = 2,
           arrow = arrow(length = unit(0.15,
                                       "inches")))

# Scan direction
CV_PtBu2Np_dir <- CVs_dir_plot$PtBu2Np +
  annotate("segment",
           x = -0.3,
           y = 2.5e-06,
           xend = 0,
           yend = 2.5e-06,
           linewidth = 2,
           arrow = arrow(length = unit(0.15,
                                       "inches")),
           color = "#00BEC4") +
  annotate("segment",
           x = 0,
           y = 4e-06,
           xend = -0.3,
           yend = 4e-06,
           linewidth = 2,
           arrow = arrow(length = unit(0.15,
                                       "inches")),
           color = "#F8766D")

# Combine with patchwork
CV_studies_PtBu2Np <- CV_PtBu2Np +
  CV_PtBu2Np_dir +
  plot_layout(axis_titles = "collect")


# PtBu3
# Regular CV
CV_PtBu3 <- CVs_pub$PtBu3 +
  annotate("segment",
           x = -0.3,
           y = 3e-06,
           xend = 0,
           yend = 3e-06,
           linewidth = 2,
           arrow = arrow(length = unit(0.15,
                                       "inches")))

# Scan direction
CV_PtBu3_dir <- CVs_dir_plot$PtBu3 +
  annotate("segment",
           x = -0.3,
           y = 3e-06,
           xend = 0,
           yend = 3e-06,
           linewidth = 2,
           arrow = arrow(length = unit(0.15,
                                       "inches")),
           color = "#00BEC4") +
  annotate("segment",
           x = 0,
           y = 4.5e-06,
           xend = -0.3,
           yend = 4.5e-06,
           linewidth = 2,
           arrow = arrow(length = unit(0.15,
                                       "inches")),
           color = "#F8766D")

# Combine with patchwork
CV_studies_PtBu3 <- CV_PtBu3 +
  CV_PtBu3_dir +
  plot_layout(axis_titles = "collect")

# Patch together scan rate plots for bidentates to make it take less space
CV_isol_dcype <- CVs_scanrate_plot_isolated_peaks$dcype +
  theme(legend.position = "bottom")

CV_isol_dppm <- CVs_scanrate_plot_isolated_peaks$dppm +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")

CV_isol_BIPHEP <- CVs_scanrate_plot_isolated_peaks$BIPHEP +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")

isol_scan_rate_bidentates <- CV_isol_dcype +
  CV_isol_dppm +
  CV_isol_BIPHEP +
  plot_layout(axis_titles = "collect")

CV_all_dppe <- CVs_scanrate_plot_all_peaks$dppe +
  theme(legend.position = "bottom")

CV_all_dcype <- CVs_scanrate_plot_all_peaks$dcype +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")

CV_all_Fdppe <- CVs_scanrate_plot_all_peaks$Fdppe +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")

full_scan_rate_bidentates <- CV_all_dppe +
  CV_all_dcype +
  CV_all_Fdppe + 
  plot_layout(axis_titles = "collect")

# Patch together the monodentate scan rate plots
mono_isol_sr <- CVs_scanrate_plot_isolated_peaks$PAd2Bu +
  (CVs_scanrate_plot_isolated_peaks$PMes3 +
     theme(axis.ticks.y = element_blank(),
           axis.text.y = element_blank(),
           legend.position = "none"))  +
  (CVs_scanrate_plot_isolated_peaks$PtBu2Cy +
     theme(axis.ticks.y = element_blank(),
           axis.text.y = element_blank())) +
  CVs_scanrate_plot_isolated_peaks$PtBu2Np +
  (CVs_scanrate_plot_isolated_peaks$PtBu3 +
     theme(axis.ticks.y = element_blank(),
           axis.text.y = element_blank())) +
  plot_layout(axis_titles = "collect",
              guides = "collect")

mono_full_sr <- CVs_scanrate_plot_all_peaks$PAd2Bu +
  (CVs_scanrate_plot_all_peaks$PtBu3+
     theme(axis.ticks.y = element_blank(),
           axis.text.y = element_blank()))  +
  plot_layout(axis_titles = "collect",
              guides = "collect")
