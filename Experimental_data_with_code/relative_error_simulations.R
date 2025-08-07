library(ggbeeswarm)

source("ligand_effect_analysis.R")

# Add weighted noise with slight variations to models
weighted_noise_fn <- function(df, no_sim){
  # Initiate seed
  seed <- 1
  
  # Initiate output df
  df_out <- df %>%
    mutate(sim_no = 0)
  
  # Initiate sequence of simulations
  sims <- seq(1, no_sim, by = 1)
  
  # Start loop of data injection
  for (n in sims) {
    # Set seed
    set.seed(seed)
    
    # Randomly sample from there
    seed_2 <- rnorm(1,0,4500)
    
    # Set that as seed now
    set.seed(seed_2)
    
    # Clone df and update logk values w/noise
    loop_df <- df %>%
      mutate(logk_rt = logk_rt 
             + rnorm(n(), 0, 0.05) * logk_rt,
             sim_no = n)
    
    # Append to output df
    df_out <- rbind(loop_df,
                    df_out)
    
    # Increase seed iteration
    seed <- seed + 1
  }
  
  return(df_out)
}

# Simulate noise injected datasets for plotting
rates_sim_error <- weighted_noise_fn(k_E_L_scaled,
                                     30) %>%
  mutate(sim_og = if_else(sim_no == 0,
                          T, F))

# Plot to visualise
rates_plot <- ggplot(data = NULL) +
  geom_beeswarm(data = rates_sim_error %>%
                  filter(sim_og == F),
                aes(y = reorder(Ligand, 
                                logk_rt),
                    x = logk_rt),
                pch = 1,
                size = 2,
                fill = "black") +
  geom_beeswarm(data = rates_sim_error %>%
                  filter(sim_og == T
                         & L_type == "bidentate"),
                aes(y = reorder(Ligand, 
                                logk_rt),
                    x = logk_rt),
                pch = 21,
                size = 4,
                fill = "#B71B1B") +
  geom_beeswarm(data = rates_sim_error %>%
                  filter(sim_og == T
                         & L_type == "monodentate"),
                aes(y = reorder(Ligand, 
                                logk_rt),
                    x = logk_rt),
                pch = 21,
                size = 4,
                fill = "#636C9D") +
  scale_y_discrete(labels = function(x) parse(
    text = rates_sim_error$Special_label_exp[
      match(x, rates_sim_error$Ligand)])) +
  scale_size_identity() +
  themething +
  labs(x = expression(paste(log[10](k),
                            " at 25 \u00B0C (DFT)")),
       y = "") +
  theme(legend.position = "none")

# Simulate for model analysis
rates_sim_error_mods <- weighted_noise_fn(k_E_L_scaled,
                                          500) %>%
  group_by(Ligand) %>%
  mutate( sim_og = if_else(sim_no == 0,
                          T, F),
         dlogk = logk_rt - logk_rt[sim_og == T]) %>%
  ungroup() %>%
  mutate(Ligand = reorder(Ligand, 
                          logk_rt))

# Alternative plot for visualising all 500 sims
sim_plot_alt <- ggplot(data = rates_sim_error_mods) +
  geom_histogram(aes(x = dlogk)) +
  facet_wrap(. ~ Ligand) +
  themething

# Function for doing model comparison
inj_nois_mod_comp_fn <- function(df){
  # Split into n dfs
  df_list <- df %>%
    group_split(sim_no)
  
  # Initialise results df
  out_df <- data.frame(model = c(),
                       sim_no = c(),
                       AICc = c(),
                       BIC = c(),
                       pA = c())
  
  # Loop through list and compare models
  for (df_sim in df_list) {
    # Create the three models
    E_model <- lm(data = df_sim, 
                       logk_rt ~ E_scaled)
    
    EV_model <- lm(data = df_sim,
                       logk_rt ~ E_scaled 
                   + Vbur_scaled)
    
    EVL_model <- lm(data = df_sim,
                            logk_rt ~ E_scaled + 
                              Vbur_scaled + L_type)
    
    # Append metrics to output df
    loop_df <- data.frame(model = c(1,2,3),
             sim_no = unique(df_sim$sim_no),
             AICc = c(AICc(E_model),
                      AICc(EV_model),
                      AICc(EVL_model)),
             BIC = c(BIC(E_model),
                      BIC(EV_model),
                      BIC(EVL_model)),
             pA = c(anova(E_model,
                          EV_model)$`Pr(>F)`[[2]],
                    anova(EV_model,
                          EVL_model)$`Pr(>F)`[[2]],
                    anova(E_model,
                          EVL_model)$`Pr(>F)`[[2]]))
    
    # Bind to output df
    out_df <- rbind(out_df,
                    loop_df)
  }
  
  return(out_df)
}

# Evaluate all simulations
mods_rates_sim_error <- inj_nois_mod_comp_fn(
  rates_sim_error_mods) %>%
  group_by(sim_no) %>%
  mutate(dAICc = AICc - min(AICc),
         dBIC = BIC - min(BIC),
         adj_dAICc = if_else(dAICc < 2,
                             0,dAICc),
         adj_dBIC = if_else(dBIC < 2,
                             0,dBIC),
         AICc_tie = if_else(length(
           unique(adj_dAICc)) < 3,
                       T, F),
         BIC_tie = if_else(length(
           unique(adj_dBIC)) < 3,
           T, F),
         pfav = if_else(pA[model == 1] < 0.05
                        & pA[model == 2] > 0.05,
                        2,
                        if_else(pA[model == 3] < 0.05
                                & pA[model == 2] < 0.05,
                                3,
                                if_else(
                                  pA[model == 1] > 0.05
                                  & pA[model == 2] > 0.05
                                  & pA[model == 3] < 0.05,
                                0,1))),
         AICcfav = if_else(AICc_tie == F
                           & dAICc == 0,
                           model,
                           if_else(AICc_tie == T 
                                   & dAICc > 2,
                                   -model,
                                   0)),
         BICfav = if_else(BIC_tie == F
                           & dBIC == 0,
                           model,
                           if_else(BIC_tie == T 
                                   & dBIC > 2,
                                   -model,
                                   0)),
         AICc_score = rank(-adj_dAICc),
         BIC_score = rank(-BIC),
         pA_score = if_else(pfav == model,
                            1,
                            if_else(pfav == 0
                                    & (model == 3
                                       | model == 1),
                                    0.5,
                                    0)),
         total_score = AICc_score + BIC_score + pA_score,
         score_rank = rank(total_score),
         total_tie = if_else(length(unique((
           score_rank))) < 3,
                             T, F),
         best_model = if_else(score_rank == 3
                              | score_rank == 2.5,
                              T,F))

# Summary table of sim results
inj_error_summary <- data.frame(
  model = c(1,2,3),
  total_best = c(nrow(mods_rates_sim_error %>%
                        filter(model == 1
                               & best_model == T)),
                 nrow(mods_rates_sim_error %>%
                        filter(model == 2
                               & best_model == T)),
                 nrow(mods_rates_sim_error %>%
                        filter(model == 3
                               & best_model == T))),
  best_no_tie = c(nrow(mods_rates_sim_error %>%
                         filter(model == 1
                                & best_model == T
                                & total_tie == F)),
                  nrow(mods_rates_sim_error %>%
                         filter(model == 2
                                & best_model == T
                                & total_tie == F)),
                  nrow(mods_rates_sim_error %>%
                         filter(model == 3
                                & best_model == T
                                & total_tie == F)))) %>%
  mutate(best_tie = total_best - best_no_tie,
         pct_best = total_best/500,
         pct_best_no_tie = best_no_tie/total_best,
         pct_best_tie = best_tie/total_best)