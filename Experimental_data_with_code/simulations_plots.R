library(ggpubr)
library(ggplot2)
library(dplyr)

source("analyse_simulations.R")
source("common_functions.R")

# Function for creating probability plots
ind_metric_make_plot <- function(df,y_aes,
                                 col_aes,
                                 y_label,
                                 leg_label) {
  
  # Get the values of ground truth model
  SD_noise_vals <- unique(df$SD_Noise)
  
  # Get the values of ground truth model
  gt_vals <- unique(df$Ground_Truth_Model)
  
  # Initialize a list of plots
  metric_plotlist <- vector("list", 
                            length = length(gt_vals))
  
  # Counter for the current position in the list
  plot_counter <- 1
  
  # Look through values of GT first
  for (gt_val in gt_vals){
    
    # Initialize plot
    metric_plot <- ggplot()
    
    # Iterate through values of SD_Noise
    for (SD_noise_val in SD_noise_vals){
      
      # Filter the dataframe for the appropriate values 
      # of SD_Noise and Ground_Truth_Model
      df_filt <- df %>%
        filter(Ground_Truth_Model == gt_val,
               SD_Noise == SD_noise_val)
      
      # Add geom_line layer
      metric_plot <- metric_plot +
        geom_line(data = df_filt,
                  aes(x = N_Obs,
                      y = {{y_aes}},
                      color = {{col_aes}}),
                  alpha = ifelse(SD_noise_val == 1, 
                                 1, 
                                 0.2),
                  linewidth = 2)
    }
    
    # Add plot title
    metric_plot <- metric_plot +
      labs(title = paste("Ground Truth Model:", gt_val),
           color = leg_label,
           y = y_label,
           x = "Sample Size") +
      ylim(0,1) +
      themething +
      scale_color_brewer(palette = "Dark2") +
      theme(plot.title = element_text(size = 10))
    
    # Add the plot to the list of plots
    metric_plotlist[[plot_counter]] <- metric_plot
    
    # Increment the counter
    plot_counter <- plot_counter + 1
  }
  
  gird <- ggarrange(metric_plotlist[[1]] +
                      theme(axis.title.x = element_text(colour = "white")),
                    metric_plotlist[[2]] +
                      theme(axis.text.y = element_blank(),
                            axis.title.y = element_text(colour = "white")),
                    metric_plotlist[[3]] +
                      theme(axis.text.y = element_blank(),
                            axis.title.y = element_text(colour = "white"),
                            axis.title.x = element_text(colour = "white")),
                    common.legend = T,
                    legend = "bottom",
                    ncol = 3,
                    nrow = 1)
  
  return(gird)
}

# Plot for how evaluating ANOVA matters
ANOVA_criteria_plot <- ind_metric_make_plot(
  ANOVA_prob,
  quotient,
  correct_col,
  "Absolute Probability Correct",
  "Selection Criterion")

# Plot for individual metrics
ind_metric_plot <- ind_metric_make_plot(
  prob_ind_metric,
  quotient,
  correct_col,
  "Absolute Probability Correct",
  "Metric")

# Filter data to look at the ones that are 
# not redundant and only when
# the metrics like AICc can actually be computed correctly
non_red_metrics <- prob_ind_metric %>%
  ungroup %>%
  filter(N_Obs > 6,
         correct_col == "AICc_correct" |
         correct_col == "BIC_correct" |
         correct_col == "ANOVA_correct")

non_red_ind_plot <- ind_metric_make_plot(
  non_red_metrics,
  quotient,
  correct_col,
  "Absolute Probability Correct",
  "Metric")

# Plot for absolute probability of being 
# correct with some metrics left out
leave_some_out_plot <- ind_metric_make_plot(
  agreement_prob,
  abs_prob_correct,
  as.character(no_par_left_out),
  "Absolute Probability Correct",
  "Number of Allowed Disagreements")

# Plot for relative probability of being 
# correct with some metrics left out
leave_some_out_plot_rel <- ind_metric_make_plot(
  agreement_prob,
  prob_correct,
  as.character(no_par_left_out),
  "Relative Probability Correct",
  "Number of Allowed Disagreements")

# Plot for probability of agreement occuring
leave_some_out_occ_plot <-ind_metric_make_plot(
  agreement_prob,
  prob_agree,
  as.character(no_par_left_out),
  "Probability of Metric Agreement",
  "Number of Allowed Disagreements")

# Plot absolute probability for non redundant metrics
non_red_abs_prob_combo_plot <- ind_metric_make_plot(
  non_red_prob,
  abs_prob_correct,
  as.character(no_par_left_out),
  "Absolute Probability Correct",
  "Number of Allowed Disagreements")

# Same plot but relative probability
non_red_rel_prob_combo_plot <- ind_metric_make_plot(
  non_red_prob,
  prob_correct,
  as.character(no_par_left_out),
  "Relative Probability Correct",
  "Number of Allowed Disagreements")

# Plot probability of non redundant metrics all agreeing
prob_agree_non_red_plot <- ind_metric_make_plot(
  non_red_prob,
  prob_agree,
  as.character(no_par_left_out),
  "Probability of Metric Agreement",
  "Number of Allowed Disagreements")

# Plot some data at three levels of noise for 
# three different ground truth
# models to show what they look like
ex_dat_1 <- generate_data(15,0.25,1,915)
ex_dat_2 <- generate_data(15,1,1,3544)
ex_dat_3 <- generate_data(15,2,1,4848)

list_ex_dfs <- list(ex_dat_1[[1]],
                    ex_dat_2[[1]],
                    ex_dat_3[[1]])

ex_dat_plot_fn <- function(list_df){
  
  list_of_plots <- list()
  
  for (df in list_df){
    
    df_plot <- ggplot(data = df) +
      geom_point(aes(x = E_ox,
                     y = logk,
                     color = L_type,
                     size = V_bur)) +
      labs(x = "Oxidation Potential (V)",
           y = expression(paste(log[10], "(k)"))) +
      scale_color_brewer(palette = "Dark2") +
      themething
    
    list_of_plots <- c(list_of_plots,
                          list(df_plot))
  }
  
  all_plots <- ggarrange(list_of_plots[[1]] +
                           theme(
                             axis.title.x = element_text(
                               colour = "white")),
                         list_of_plots[[2]] +
                           theme(axis.text.y = element_blank(),
                                 axis.title.y = element_text(
                                   colour = "white")),
                         list_of_plots[[3]] +
                           theme(axis.text.y = element_blank(),
                                 axis.title.y = element_text(
                                   colour = "white"),
                                 axis.title.x = element_text(
                                   colour = "white")),
                         nrow = 1,
                         ncol = 3,
                         common.legend = T,
                         legend = "bottom")
  
  return(all_plots)
}

ex_sim_data_plots <- ex_dat_plot_fn(list_ex_dfs)

# Plot for probability correct based on observed consensus
# Function for creating probability plots
consensus_make_plot <- function(df) {
  
  # Get the values of ground truth model
  SD_noise_vals <- unique(df$SD_Noise)
  
  # Get the values of ground truth model
  gt_vals <- unique(df$no_par_left_out)
  
  # Initialize a list of plots
  metric_plotlist <- vector("list", 
                            length = length(gt_vals))
  
  # Counter for the current position in the list
  plot_counter <- 1
  
  # Look through values of GT first
  for (gt_val in gt_vals){
    
    # Initialize plot
    metric_plot <- ggplot()
    
    # Iterate through values of SD_Noise
    for (SD_noise_val in SD_noise_vals){
      
      # Filter the dataframe for the appropriate values 
      # of SD_Noise and Ground_Truth_Model
      df_filt <- df %>%
        filter(no_par_left_out == gt_val,
               SD_Noise == SD_noise_val)
      
      # Add geom_line layer
      metric_plot <- metric_plot +
        geom_line(data = df_filt,
                  aes(x = N_Obs,
                      y = prob_correct,
                      color = as.character(metric_consensus)),
                  alpha = ifelse(SD_noise_val == 1, 
                                 1, 
                                 0.2),
                  linewidth = 2)
    }
    
    # Add plot title
    metric_plot <- metric_plot +
      labs(title = paste("Number of Disagreements:", gt_val),
           color = "Metric Consensus",
           y = "Probability Correct",
           x = "Sample Size") +
      ylim(0,1) +
      themething +
      scale_color_manual(values = c("#636C9D",
                                     "#B71B1B",
                                     "#4E525B")) +
      theme(plot.title = element_text(size = 10))
    
    # Add the plot to the list of plots
    metric_plotlist[[plot_counter]] <- metric_plot
    
    # Increment the counter
    plot_counter <- plot_counter + 1
  }
  
  gird <- ggarrange(metric_plotlist[[1]],
                    metric_plotlist[[2]] +
                      theme(axis.text.y = element_blank(),
                            axis.title.y = element_text(
                              colour = "white")),
                    common.legend = T,
                    legend = "bottom",
                    ncol = 2,
                    nrow = 1)
  
  return(gird)
}

consensus_probability_plot <- consensus_make_plot(
  consensus_prob)
