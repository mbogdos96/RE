library(dplyr)
library(ggplot2)
library(ggrepel)
library(AICcmodavg)
library(patchwork)
library(lavaan)

source("Eox_descriptors.R")
source("literature_data.R")

#==============
# Data handling
#==============
# Palladacycles
# Merge ligand HOMO energies with df containing 
# other info
all_data_df <- merge(k_E_L_scaled,
                     ligand_HOMO_df) %>%
  left_join(pdcy_free_L_df %>%
              dplyr::select(c("Ligand",
                              "Potential_free_L")),
            by = "Ligand") %>%
  mutate(L_HOMO_scaled = HOMO_energy_kcal_mol,
         scaled_logk = logk_rt,
         scaled_free_L_E = Potential_free_L) %>%
  mutate(across(c(L_HOMO_scaled, 
                scaled_logk,
                scaled_free_L_E),
              ~ scale(.,
                      scale = T,
                      center = T)))

# Raw estimates palladacycle data
pdcy_df_raw <- all_data_df %>%
  rename(E = "Potential",
         CN = "L_type",
         P_don = "Potential_free_L") %>%
  mutate(CN = relevel(CN,
                      ref = "bidentate"))

# Scaled estimates palladacycle data
pdcy_df_scaled <- all_data_df %>%
  dplyr::select(-c(Vbur,logk_rt)) %>%
  rename(E = "E_scaled",
         CN = "L_type",
         P_don = "scaled_free_L_E",
         Vbur = "Vbur_scaled",
         logk_rt = "scaled_logk") %>%
  mutate(CN = relevel(CN,
                      ref = "bidentate"))

# Literature
# Raw estimates literature complexes
lit_raw_df <- read.csv(file.path(
  "./literature_data/inference_subsample_lit.csv"
)) %>%
  rename(E = "Soln_VIP_eV",
         Vbur = "V_bur",
         CN = "Coord_no",
         P_don = "Ligand.Ehomo...Eh",
         Nu_par = "Nucleophilicity_eV",
         E_par = "Electrophilicity_eV",
         E_don = "Electrophile.Anion.HOMO...eV") %>%
  mutate(CN = as.factor(CN),
         CN = relevel(CN,
                      ref = "4")) %>%
  left_join(subsample_lit_df %>%
              dplyr::select(c(ID,
                              logk_rt)),
            by = "ID")

# Scaled estimates literature complexes
lit_scaled_df <- lit_raw_df %>%
  mutate(across(c(E, Vbur, P_don, Nu_par,
                  E_par, E_don, logk_rt),
                ~ scale(.,
                      scale = T,
                      center = T)))

#================================
# Mediation analysis calculations
#================================
# Function for applying models and extracting 
# path coefficients
med_fn <- function(inp_df, 
                   app_mod) {
  
  # Apply the linear models based on app_mod
  if (app_mod == "partial") {
    # Apply k pred lm()
    lm_E_S_logk <- lm(logk_rt ~ E + Vbur, 
                      data = inp_df)
    
    # Apply mediation lm()
    lm_CN_S_E <- lm(E ~ CN + Vbur + P_don, 
                    data = inp_df)
  } else {
    # Apply k pred lm()
    lm_E_S_logk <- lm(logk_rt ~ E + Vbur 
                      + Nu_par + E_par, 
                      data = inp_df)
    
    # Apply mediation lm()
    lm_CN_S_E <- lm(E ~ CN + Vbur + P_don 
                    + Nu_par + E_don, 
                    data = inp_df)
  }
  
  # Construct data frame with estimates and errors
  est_df <- data.frame(summary(lm_E_S_logk)$coefficients) %>%
    head(3) %>%
    rbind(data.frame(summary(lm_CN_S_E)$coefficients) %>%
            head(4)) %>%
    rename(Std_error = "Std..Error") %>%
    dplyr::select(c(Estimate, Std_error)) %>%
    mutate(Est_upr = Estimate + Std_error,
           Est_lwr = Estimate - Std_error,
           effect = c("", "e", "s", "", 
                      "s^cn", "s^i", "d")) %>%
    filter(effect != "")
  
  # Define effects and operations
  effects <- c("s^e_k", "s^cn_k", 
               "s^tot", "s^tot_CN4")
  operations <- list(
    s_e_k = function(e, 
                     s_i) e * s_i,
    s_cn_k = function(e, 
                      s_cn) e * s_cn,
    s_tot = function(e, s, 
                     s_i, s_cn) 
      (e * s_i) + (e * s_cn) + s,
    s_tot_CN4 = function(e, 
                         s, s_i) (e * s_i) + s
  )
  
  # Standard error propagation functions for each effect
  error_propagation <- list(
    s_e_k = function(
    e, s_i, 
    se_e, 
    se_s_i) sqrt((s_i * se_e)^2 + (e * se_s_i)^2),
    s_cn_k = function(
    e, s_cn, 
    se_e, se_s_cn) sqrt((s_cn * se_e)^2 + (e * se_s_cn)^2),
    s_tot = function(
    e, s, s_i, s_cn, 
    se_e, se_s, 
    se_s_i, 
    se_s_cn) 
      sqrt((s_i * se_e)^2 + (e * se_s_i)^2 + (s_cn * se_e)^2 + (e * se_s_cn)^2 + se_s^2),
    s_tot_CN4 = function(e, s, s_i, 
                         se_e, se_s, se_s_i) 
      sqrt((s_i * se_e)^2 + (e * se_s_i)^2 + se_s^2)
  )
  
  # Initialize lists for new rows and standard errors
  new_rows <- list()
  new_errors <- list()
  
  # Loop over Estimate to compute values for each effect
  for (col in c("Estimate", "Std_error")) {
    e <- est_df[est_df$effect == "e", "Estimate"]
    s <- est_df[est_df$effect == "s", "Estimate"]
    s_i <- est_df[est_df$effect == "s^i", "Estimate"]
    s_cn <- est_df[est_df$effect == "s^cn", "Estimate"]
    
    se_e <- est_df[est_df$effect == "e", "Std_error"]
    se_s <- est_df[est_df$effect == "s", "Std_error"]
    se_s_i <- est_df[est_df$effect == "s^i", "Std_error"]
    se_s_cn <- est_df[est_df$effect == "s^cn", "Std_error"]
    
    # Calculate new values for each effect based on Estimate
    if (col == "Estimate") {
      new_vals <- c(
        operations$s_e_k(e, s_i),
        operations$s_cn_k(e, s_cn),
        operations$s_tot(e, s, s_i, s_cn),
        operations$s_tot_CN4(e, s, s_i)
      )
      new_rows[[col]] <- new_vals
    }
    
    # Calculate propagated standard errors based on Std_error
    if (col == "Std_error") {
      new_se_vals <- c(
        error_propagation$s_e_k(
          e, s_i, se_e, se_s_i),
        error_propagation$s_cn_k(
          e, s_cn, se_e, se_s_cn),
        error_propagation$s_tot(
          e, s, s_i, s_cn, se_e, se_s, se_s_i, se_s_cn),
        error_propagation$s_tot_CN4(
          e, s, s_i, se_e, se_s, se_s_i)
      )
      new_errors[[col]] <- new_se_vals
    }
  }
  
  # Create the new rows data frame
  full_est_df <- data.frame(
    Estimate = new_rows$Estimate,
    Std_error = new_errors$Std_error,
    Est_upr = new_rows$Estimate + new_errors$Std_error,
    Est_lwr = new_rows$Estimate - new_errors$Std_error,
    effect = effects
  )
  
  # Combine the data frames
  all_est_df <- rbind(est_df, full_est_df) %>%
    mutate(type_effect = ifelse(effect %in% c(
      "e", "s^tot", "s^tot_CN4"),
                                "Total Effects",
                                ifelse(effect %in% c("s^cn_k", "s^e_k", "s"),
                                       "Effects on Rate",
                                       "Indirect Effects")))
  
  # Return the final data frame with estimates
  return(all_est_df)
}

# Apply mediation analysis function 
# Raw estimates of palladacycle data
med_raw_pdcy_df <- med_fn(pdcy_df_raw,
                          "partial") %>%
  mutate(scaling = "N",
         sample = "PdCy")

# Scaled palladacycle data
med_scaled_pdcy_df <- med_fn(pdcy_df_scaled,
                             "partial")%>%
  mutate(scaling = "Y",
         sample = "PdCy")

# Raw lit subsample values
med_raw_lit_df <- med_fn(lit_raw_df,
                         "full") %>%
  mutate(scaling = "N",
         sample = "Lit")

# Scaled lit subsampled values
med_scaled_lit_df <- med_fn(lit_scaled_df,
                            "full") %>%
  mutate(scaling = "Y",
         sample = "Lit")

# Merge all into one df
med_master_df <- rbind(med_scaled_pdcy_df,
                       med_raw_lit_df,
                       med_scaled_lit_df,
                       med_raw_pdcy_df) %>%
  filter(effect != "d") %>%
  group_by(scaling,
           sample) %>%
  mutate(est_rel_e = Estimate/Estimate[effect == "e"],
         est_upr_rel_e = Est_upr/Estimate[effect == "e"],
         est_lwr_rel_e = Est_lwr/Estimate[effect == "e"])

#=========
# Plotting
#=========
# Function for making the estimate plots
est_plot_fn <- function(in_df,
                        aes_y,
                        error_y_max,
                        error_y_min,
                        scaling_vec,
                        sample_vec,
                        y_lab_lab,
                        color_by){
  
  # Filter input df
  plot_df <- in_df %>%
    filter(scaling %in% scaling_vec) %>%
    filter(sample %in% sample_vec)
  
  # Make the effect names expressions
  med_an_x_vec <- c(expression(e),
                    expression(s^{d}),
                    expression(s^{cn}),
                    expression({s}^{i}),
                    expression(s[k]^{e}),
                    expression(s[k]^{cn}),
                    expression(s[CN3]^{tot}),
                    expression(s[CN4]^{tot}))
  
  names(med_an_x_vec) <- unique(med_master_df$effect)
  
  # Make the plot
  out_plot <- ggplot(data = plot_df,
         aes(x = factor(effect,
                        levels = c("e",
                                   "s^tot",
                                   "s^tot_CN4",
                                   "s^cn_k",
                                   "s^e_k",
                                   "s",
                                   "s^cn",
                                   "s^i",
                                   "d")),
             fill = {{color_by}})) +
    geom_hline(yintercept = 0, 
               linetype = "dashed", 
               color = "#4E525A", 
               linewidth = 1,
               alpha = 0.3) +
    geom_errorbar(aes(ymin = {{error_y_min}},
                      ymax = {{error_y_max}},
                      color = {{color_by}}),
                  width = 0,
                  linewidth = 1.25,
                  alpha = 0.3,
                  position = position_dodge(0.5)) +
    geom_point(aes(y = {{aes_y}}),
               position = position_dodge(0.5),
               size = 3,
               pch = 21) +
    facet_grid(cols = vars(factor(type_effect,
                                  levels = c(
                                    "Total Effects",
                                    "Effects on Rate",
                                    "Indirect Effects"))), 
               scales = "free") +
    labs(x = "Effect",
         y = y_lab_lab) +
    scale_fill_manual(values = c("#B71B1B",
                                 "#636C9D")) +
    scale_color_manual(values = c("#B71B1B",
                                  "#636C9D")) +
    scale_x_discrete(labels = med_an_x_vec) +
    themething +
    theme(legend.position = "none")
  
  return(out_plot)
}

# Only scaled pdcy values relative to e
est_sc_pdcy_rel_plot <- est_plot_fn(
  med_master_df,
  est_rel_e,
  est_upr_rel_e,
  est_lwr_rel_e,
  c("Y"),
  c("PdCy"),
  "Estimate Relative to e",
  sample)

# Plot of raw vs scaled estimates for palladacycles
# relative to e
est_rel_raw_v_sc_pdcy_plot <- est_plot_fn(
  med_master_df %>%
    mutate(scaling = as.factor(scaling),
           scaling = relevel(scaling,
                            ref = "Y")),
  est_rel_e,
  est_upr_rel_e,
  est_lwr_rel_e,
  c("Y","N"),
  c("PdCy"),
  "Estimate Relative to e",
  scaling) +
  theme(legend.position = "bottom") +
  labs(fill = "Data scaling",
       color = "Data scaling")

# Plot of raw vs scaled for palladacycles in
# absolute numbers
est_raw_v_sc_pdcy_plot <- est_plot_fn(
  med_master_df %>%
    mutate(scaling = as.factor(scaling),
           scaling = relevel(scaling,
                             ref = "Y")),
  Estimate,
  Est_upr,
  Est_lwr,
  c("Y","N"),
  c("PdCy"),
  "Estimate",
  scaling) +
  theme(legend.position = "bottom") +
  labs(fill = "Data scaling",
       color = "Data scaling")

# Combine the raw vs scaled plots
sc_v_raw_comb_plot <- est_rel_raw_v_sc_pdcy_plot +
  est_raw_v_sc_pdcy_plot +
  plot_layout(axes = "collect",
              guides = "collect") &
  theme(legend.position = "bottom")

# Lit estimates path coef raw v scaled
lit_raw_v_scaled <- est_plot_fn(
  med_master_df %>%
    mutate(scaling = as.factor(scaling),
           scaling = relevel(scaling,
                             ref = "Y")),
  Estimate,
  Est_upr,
  Est_lwr,
  c("Y","N"),
  c("Lit"),
  "Estimate",
  scaling) +
  theme(legend.position = "bottom") +
  labs(fill = "Data scaling",
       color = "Data scaling")

# Lit estimates path coef raw v scaled rel values
lit_raw_v_scaled_rel <- est_plot_fn(
  med_master_df %>%
    mutate(scaling = as.factor(scaling),
           scaling = relevel(scaling,
                             ref = "Y")),
  est_rel_e,
  est_upr_rel_e,
  est_lwr_rel_e,
  c("Y","N"),
  c("Lit"),
  "Estimate",
  scaling) +
  theme(legend.position = "bottom") +
  labs(fill = "Data scaling",
       color = "Data scaling")

# Combine lit plots
lit_raw_v_scal_comb <- lit_raw_v_scaled_rel +
  lit_raw_v_scaled +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Plot of both lit and pdcy scaled values, relative to e
est_both_rel_plot <- est_plot_fn(
  med_master_df %>%
    mutate(sample = as.factor(sample),
           sample = relevel(sample,
                            ref = "PdCy")),
  est_rel_e,
  est_upr_rel_e,
  est_lwr_rel_e,
  c("Y"),
  c("Lit",
    "PdCy"),
  "Estimate Relative to e",
  sample) +
  theme(legend.position = "bottom") +
  labs(fill = "Dataset",
       color = "Dataset")

# Plot of both lit and pdcy scaled values
est_both_plot <- est_plot_fn(
  med_master_df %>%
    mutate(sample = as.factor(sample),
           sample = relevel(sample,
                            ref = "PdCy")),
  Estimate,
  Est_upr,
  Est_lwr,
  c("Y"),
  c("Lit",
    "PdCy"),
  "Estimate relative to e",
  sample) +
  theme(legend.position = "bottom") +
  labs(fill = "Dataset",
       color = "Dataset")

# Combine plots comparing scaled lit and scaled pdcy
lit_v_pdcy_plot <- est_both_rel_plot + 
  est_both_plot +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

#=====================================
# Validation w/SEM from lavaan package
#=====================================
# Define the limited DAG
DAG_structure_partial <- '
# Direct effects on rate
logk_rt ~ E + Vbur

# Effect on E
E ~ Vbur + CN + P_don
'

# Use SEM on limited DAG
sem_partial_pdcy <- sem(DAG_structure_partial,
                        data = pdcy_df_scaled)

# Define the limited DAG w/composite effects
DAG_structure_partial_comp_eff <- '
# Direct effects on rate
logk_rt ~ e*E + sd*Vbur

# Effect on E
E ~ si*Vbur + scn*CN + d*P_don

# Indirect effect S on k
sek := si*e

# Effect of coord no change on k
scnk := scn*e

# Total effect sterics T-shaped
stotCN3 := (si*e) + sd

# Total effect sterics square-planar
stotCN4 := (scn*e) + (si*e) + sd
'

# Use SEM on limited DAG
sem_partial_pdcy_err <- sem(DAG_structure_partial_comp_eff,
                        data = pdcy_df_scaled,
                        se = "boot",
                        bootstrap = 500)

# Define the expanded DAG for lit
DAG_structure_full <- '
logk_rt ~ E + Vbur + Nu_par + E_par
E ~ Vbur + CN + P_don + Nu_par + E_don
'

# Use SEM on lit data
sem_full_lit <- sem(DAG_structure_full,
                    data = lit_scaled_df)

#==========================================
# Reviewer 2 fractional coordination number
#==========================================
# Function for applying models and extracting 
# path coefficients
med_fn_frCN <- function(inp_df, 
                   app_mod) {
  
  # Apply the linear models based on app_mod
  if (app_mod == "partial") {
    # Apply k pred lm()
    lm_E_S_logk <- lm(logk_rt ~ E + Vbur, 
                      data = inp_df)
    
    # Apply mediation lm()
    lm_CN_S_E <- lm(E ~ CN + Vbur + P_don, 
                    data = inp_df)
  } else {
    # Apply k pred lm()
    lm_E_S_logk <- lm(logk_rt ~ E + Vbur 
                      + Nu_par + E_par, 
                      data = inp_df)
    
    # Apply mediation lm()
    lm_CN_S_E <- lm(E ~ CN + Vbur + P_don 
                    + Nu_par + E_don, 
                    data = inp_df)
  }
  
  # Construct data frame with estimates and errors
  est_df <- data.frame(summary(lm_E_S_logk)$coefficients) %>%
    head(3) %>%
    rbind(data.frame(summary(lm_CN_S_E)$coefficients) %>%
            head(4)) %>%
    rename(Std_error = "Std..Error") %>%
    dplyr::select(c(Estimate, Std_error)) %>%
    mutate(Est_upr = Estimate + Std_error,
           Est_lwr = Estimate - Std_error,
           effect = c("", "e", "s", "", 
                      "s^cn", "s^i", "d")) %>%
    filter(effect != "")
  
  # Define effects and operations
  effects <- c("s^e_k", "s^cn_k", 
               "s^tot", "s^tot_CN4")
  operations <- list(
    s_e_k = function(e, 
                     s_i) e * s_i,
    s_cn_k = function(e, 
                      s_cn) e * s_cn,
    s_tot = function(e, s, 
                     s_i, s_cn) 
      (e * s_i) - (e * s_cn) + s,
    s_tot_CN4 = function(e, 
                         s, s_i) (e * s_i) + s
  )
  
  # Standard error propagation functions for each effect
  error_propagation <- list(
    s_e_k = function(
    e, s_i, 
    se_e, 
    se_s_i) sqrt((s_i * se_e)^2 + (e * se_s_i)^2),
    s_cn_k = function(
    e, s_cn, 
    se_e, se_s_cn) sqrt((s_cn * se_e)^2 + (e * se_s_cn)^2),
    s_tot = function(
    e, s, s_i, s_cn, 
    se_e, se_s, 
    se_s_i, 
    se_s_cn) 
      sqrt((s_i * se_e)^2 + (e * se_s_i)^2 + (s_cn * se_e)^2 + (e * se_s_cn)^2 + se_s^2),
    s_tot_CN4 = function(e, s, s_i, 
                         se_e, se_s, se_s_i) 
      sqrt((s_i * se_e)^2 + (e * se_s_i)^2 + se_s^2)
  )
  
  # Initialize lists for new rows and standard errors
  new_rows <- list()
  new_errors <- list()
  
  # Loop over Estimate to compute values for each effect
  for (col in c("Estimate", "Std_error")) {
    e <- est_df[est_df$effect == "e", "Estimate"]
    s <- est_df[est_df$effect == "s", "Estimate"]
    s_i <- est_df[est_df$effect == "s^i", "Estimate"]
    s_cn <- est_df[est_df$effect == "s^cn", "Estimate"]
    
    se_e <- est_df[est_df$effect == "e", "Std_error"]
    se_s <- est_df[est_df$effect == "s", "Std_error"]
    se_s_i <- est_df[est_df$effect == "s^i", "Std_error"]
    se_s_cn <- est_df[est_df$effect == "s^cn", "Std_error"]
    
    # Calculate new values for each effect based on Estimate
    if (col == "Estimate") {
      new_vals <- c(
        operations$s_e_k(e, s_i),
        operations$s_cn_k(e, s_cn),
        operations$s_tot(e, s, s_i, s_cn),
        operations$s_tot_CN4(e, s, s_i)
      )
      new_rows[[col]] <- new_vals
    }
    
    # Calculate propagated standard errors based on Std_error
    if (col == "Std_error") {
      new_se_vals <- c(
        error_propagation$s_e_k(
          e, s_i, se_e, se_s_i),
        error_propagation$s_cn_k(
          e, s_cn, se_e, se_s_cn),
        error_propagation$s_tot(
          e, s, s_i, s_cn, se_e, se_s, se_s_i, se_s_cn),
        error_propagation$s_tot_CN4(
          e, s, s_i, se_e, se_s, se_s_i)
      )
      new_errors[[col]] <- new_se_vals
    }
  }
  
  # Create the new rows data frame
  full_est_df <- data.frame(
    Estimate = new_rows$Estimate,
    Std_error = new_errors$Std_error,
    Est_upr = new_rows$Estimate + new_errors$Std_error,
    Est_lwr = new_rows$Estimate - new_errors$Std_error,
    effect = effects
  )
  
  # Combine the data frames
  all_est_df <- rbind(est_df, 
                      full_est_df) %>%
    mutate(type_effect = ifelse(effect %in% c(
      "e", "s^tot", "s^tot_CN4"),
      "Total Effects",
      ifelse(effect %in% c("s^cn_k", "s^e_k", "s"),
             "Effects on Rate",
             "Indirect Effects")))
  
  # Return the final data frame with estimates
  return(all_est_df)
}

# Merge the data for scaled data
med_fr_CN <- read.csv(
  "./complex_descriptors/fractional_CN.csv") %>%
  merge(pdcy_df_scaled) %>%
  mutate(CN = scale(Pd_fractional_CN,
                                  scale = T,
                                  center = T))

# Define the DAG
DAG_fr_CN <- '
# Direct effects on rate
logk_rt ~ E + Vbur

# Effect on E
E ~ Vbur + CN + P_don
'

# Use SEM on limited DAG
sem_fr_CN <- sem(DAG_fr_CN,
                 data = med_fr_CN)

# Merge the data for raw data
med_fr_CN_raw <- read.csv(
  "./complex_descriptors/fractional_CN.csv") %>%
  merge(pdcy_df_raw) %>%
  mutate(CN = Pd_fractional_CN)

# Use mediation function on fractional coordination 
# number data
med_est_frac_CN <- med_fn_frCN(med_fr_CN,
                          "partial") %>%
  mutate(scaling = "Y") %>%
  rbind(med_fn_frCN(med_fr_CN_raw,
               "partial") %>%
          mutate(scaling = "N")) %>%
  filter(effect != "d") %>%
  group_by(scaling) %>%
  mutate(est_rel_e = Estimate/Estimate[effect == "e"],
         est_upr_rel_e = Est_upr/Estimate[effect == "e"],
         est_lwr_rel_e = Est_lwr/Estimate[effect == "e"],
         sample = "PdCy")

# Plot scaled vs not scaled path coefs rel to e
est_frac_CN_plot <- est_plot_fn(
  med_est_frac_CN %>%
    mutate(scaling = as.factor(scaling),
           scaling = relevel(scaling,
                             ref = "Y")),
  est_rel_e,
  est_upr_rel_e,
  est_lwr_rel_e,
  c("Y","N"),
  c("PdCy"),
  "Estimate Relative to e",
  scaling) +
  theme(legend.position = "bottom") +
  labs(fill = "Data scaling",
       color = "Data scaling")

# Absolute numbers for fractional coordination number
est_raw_v_sc_pdcy_plot_frCN <- est_plot_fn(
  med_est_frac_CN %>%
    mutate(scaling = as.factor(scaling),
           scaling = relevel(scaling,
                             ref = "Y")),
  Estimate,
  Est_upr,
  Est_lwr,
  c("Y","N"),
  c("PdCy"),
  "Estimate",
  scaling) +
  theme(legend.position = "bottom") +
  labs(fill = "Data scaling",
       color = "Data scaling")

# Fractional CN combined plot
frCN_comb_plot <- est_frac_CN_plot +
  est_raw_v_sc_pdcy_plot_frCN

#==============
# MS and Review
#==============
# Scaled only for Fig 10
est_frac_CN_plot_scaled <- est_plot_fn(
  med_est_frac_CN %>%
    filter(scaling == "Y"),
  est_rel_e,
  est_upr_rel_e,
  est_lwr_rel_e,
  c("Y"),
  c("PdCy"),
  "Estimate Relative to e",
  scaling)

# Combine the plots for scaled relative estimates
# for categorical and continuous CN
rev2_plots_estimates <- est_sc_pdcy_rel_plot +
  est_frac_CN_plot_scaled
