library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)

source("common_functions.R")

# Import data
xrd_structural_data <- read.csv(
  file.path("./complex_descriptors/Pd_Xray.csv"))

# Bond lengths
xrd_lengths <- xrd_structural_data %>%
  mutate(Code = as.numeric(Code)) %>%
  dplyr::select(Name,
                Code,
                Structure_source,
                Coordination_Number,
                contains("_length"))

# Plot for bond lengths
xrd_PdC_PdN_l <- ggplot(data = xrd_lengths) +
  geom_label_repel(data = . %>%
                     filter(Structure_source == "New",
                            Pd.C_length < 2.025),
                   aes(x = Pd.N_length,
                       y = Pd.C_length,
                       label = Code),
                   nudge_x = 0.005,
                   nudge_y = -0.008) +
  geom_label_repel(data = . %>%
                     filter(Structure_source == "New",
                            Pd.C_length > 2.025),
                   aes(x = Pd.N_length,
                       y = Pd.C_length,
                       label = Code),
                   nudge_y = 0.007) +
  geom_point(aes(x = Pd.N_length,
                 y = Pd.C_length,
                 color = as.character(Coordination_Number)),
             size = 4) +
  scale_color_manual(values = c("#636C9D",
                                "#B71B1B")) +
  labs(x = expression(paste("Pd-N Bond Length (",
                            ring(A),
                            ")")),
       y = expression(paste("Pd-C Bond Length (",
                                ring(A),
                                ")"))) +
  themething +
  theme(legend.position = "none")

# Second plot of bond lengths
xrd_PdP_CN_l <- ggplot(data = xrd_lengths) +
  geom_label_repel(data = . %>%
                     filter(Coordination_Number == 3 
                            & Structure_source == "New"),
                   aes(x = Pd.P_transN_length,
                       y = C..N_length,
                       label = Code),
                   nudge_y = 0.02) +
  geom_label_repel(data = . %>%
                     filter(
                       Coordination_Number == 4
                       & Structure_source == "New"
                       & Name %in% c(
                         "F[4]*-N^{Cl}*Ms-OMe-Fdppe",
                         "F[4]*-N^{Cl}*Ms-OMe-dcype")),
                   aes(x = Pd.P_transN_length,
                       y = C..N_length,
                       label = Code),
                   nudge_y = 0.035,
                   nudge_x = -0.01) +
  geom_label_repel(data = . %>%
                     filter(
                       Coordination_Number == 4
                       & Structure_source == "New"
                       & !(Name %in% c(
                         "F[4]*-N^{Cl}*Ms-OMe-Fdppe",
                         "F[4]*-N^{Cl}*Ms-OMe-dcype"))),
                   aes(x = Pd.P_transN_length,
                       y = C..N_length,
                       label = Code),
                   nudge_y = -0.025,
                   nudge_x = 0.01) +
  geom_point(aes(x = Pd.P_transN_length,
                 y = C..N_length,
                 color = as.character(Coordination_Number)),
             size = 4) +
  scale_color_manual(values = c("#636C9D",
                                "#B71B1B")) +
  labs(x = expression(paste("Pd-",
                            {P}^{trans-N},
                            " Bond Length (",
                            ring(A),
                            ")")),
       y = expression(paste("C-N Distance (",
                            ring(A),
                            ")"))) +
  themething +
  theme(legend.position = "none")

# Combine plots of bond lengths
xrd_lengths_plot <- xrd_PdC_PdN_l +
  xrd_PdP_CN_l

# Bond angles
# Simple bond angle analysis
xrd_angles <- xrd_structural_data %>%
  dplyr::select(Name,
                Code,
                Structure_source,
                Coordination_Number,
                contains("_angle"))

# First bond angle plot
xrd_PtNPdC_CPdN_a <- ggplot(data = xrd_angles) +
  geom_label_repel(data = . %>%
                     filter(Structure_source == "New"
                            & C.Pd.N_angle < 84),
                   aes(x = P_transN.Pd.C_angle,
                       y = C.Pd.N_angle,
                       label = Code),
                   nudge_x = 2,
                   nudge_y = -1) +
  geom_label_repel(data = . %>%
                     filter(Structure_source == "New"
                            & C.Pd.N_angle > 84),
                   aes(x = P_transN.Pd.C_angle,
                       y = C.Pd.N_angle,
                       label = Code),
                   nudge_x = 1.5,
                   nudge_y = 1) +
  geom_point(aes(x = P_transN.Pd.C_angle,
                 y = C.Pd.N_angle,
                 color = as.character(Coordination_Number)),
             size = 4) +
  scale_color_manual(values = c("#636C9D",
                                "#B71B1B")) +
  labs(x = expression(paste("C-Pd-",
                            {P}^{trans-N},
                            " Bond Angle (\u00B0)")),
       y = "C-Pd-N Bond Angle (\u00B0)") +
  themething +
  theme(legend.position = "none")

# Second angle plot
xrd_PtNPdPtC_PtCPdN_a <- ggplot(data = xrd_angles) +
  geom_label_repel(data = . %>%
                     filter(Structure_source == "New"
                            & P_transN.Pd.P_transC_angle > 95),
                   aes(x = P_transN.Pd.P_transC_angle,
                       y = P_transC.Pd.N_angle,
                       label = Code),
                   nudge_x = -1.25,
                   nudge_y = -1) +
  geom_label_repel(data = . %>%
                     filter(Structure_source == "New"
                            & P_transN.Pd.P_transC_angle < 95),
                   aes(x = P_transN.Pd.P_transC_angle,
                       y = P_transC.Pd.N_angle,
                       label = Code),
                   nudge_x = 2,
                   nudge_y = 1) +
  geom_point(aes(x = P_transN.Pd.P_transC_angle,
                 y = P_transC.Pd.N_angle,
                 color = as.character(Coordination_Number)),
             size = 4) +
  scale_color_manual(values = c("#636C9D",
                                "#B71B1B")) +
  labs(y = expression(paste("N-Pd-",
                            {P}^{trans-C},
                            " Bond Angle (\u00B0)")),
       x = expression(paste({P}^{trans-N},
                            "-Pd-",
                            {P}^{trans-C},
                            " Bond Angle (\u00B0)"))) +
  themething +
  theme(legend.position = "none")

# Combine angles plots
xrd_angles_plot <- xrd_PtNPdC_CPdN_a +
  xrd_PtNPdPtC_PtCPdN_a

# Deviation from planarity
# t4 and t4' for square planar
xrd_t4_df <- xrd_angles %>%
  filter(Coordination_Number == 4) %>%
  dplyr::select(-c(CPdP.CPdN_plane_angle,
                Pd.P_transN.C1_angle,
                Pd.P_transN.C2_angle,
                Pd.P_transN.C3.bridge._angle,
                C1.P_transN.C2_angle,
                C1.P_transN.C3.bridge._angle,
                C2.P_transN.C3.bridge._angle,
                Pd.P_transC.C1_angle,
                Pd.P_transC.C2_angle,
                Pd.P_transC.C3.bridge._angle,
                C1.P_transC.C2_angle,
                C1.P_transC.C3.bridge._angle,
                C2.P_transC.C3.bridge._angle)) %>%
  mutate(Code = as.character(Code)) %>%
  rowwise() %>%
  mutate(alpha_angle = sort(c(P_transN.Pd.C_angle,
                        P_transN.Pd.N_angle,
                        P_transN.Pd.P_transC_angle,
                        P_transC.Pd.C_angle,
                        P_transC.Pd.N_angle,
                        C.Pd.N_angle),
                      decreasing = T)[1],
         beta_angle = sort(c(P_transN.Pd.C_angle,
                      P_transN.Pd.N_angle,
                      P_transN.Pd.P_transC_angle,
                      P_transC.Pd.C_angle,
                      P_transC.Pd.N_angle,
                      C.Pd.N_angle),
                      decreasing = T)[2]) %>%
  ungroup() %>%
  mutate(tau_4 = (360 - (alpha_angle + beta_angle))/141,
         tau_4_prime = (beta_angle - alpha_angle)/(360 - 109.5) 
         + (180 - beta_angle)/(180 - 109.5))

# tau4 plot
xrd_tau_plot <- ggplot(data = xrd_t4_df) +
  geom_abline(slope = 1,
              intercept = 0,
              color = "#4E525A",
              linewidth = 1) +
  geom_label_repel(data = . %>%
                     filter(Structure_source == "New"
                            & tau_4 > 0.2),
                   aes(x = tau_4,
                       y = tau_4_prime,
                       label = Code),
                   nudge_x = -0.02,
                   nudge_y = 0.01) +
  geom_label_repel(data = . %>%
                     filter(Structure_source == "New"
                            & tau_4 < 0.2),
                   aes(x = tau_4,
                       y = tau_4_prime,
                       label = Code),
                   nudge_x = 0.02,
                   nudge_y = -0.01) +
  geom_point(aes(x = tau_4,
                 y = tau_4_prime,
                 color = as.character(Coordination_Number)),
             size = 4) +
  scale_color_manual(values = c("#B71B1B",
                                "#636C9D")) +
  labs(x = expression(paste(tau[4],
                      " Geometry Index")),
       y = expression(paste(tau[4],
                            "' Geometry Index"))) +
  themething +
  theme(legend.position = "none")

# Plane deviation for T-shaped
xrd_plane_angle <- ggplot(data = xrd_angles) +
  geom_label_repel(data = . %>%
                     filter(Structure_source == "New"),
                   aes(x = Structure_source,
                       y = CPdP.CPdN_plane_angle,
                       label = Code),
                   nudge_y = 1) +
  geom_point(aes(x = Structure_source,
                 y = CPdP.CPdN_plane_angle,
                 color = as.character(Coordination_Number)),
             size = 4) +
  scale_color_manual(values = c("#636C9D",
                                "#B71B1B")) +
  labs(y = expression(paste("Angle between the N-Pd-C and C-Pd-",
                            {P}^{trans-N},
                            " Planes (\u00B0)")),
       x = "Structure") +
  themething +
  theme(legend.position = "none")

# combine tau4 and deviation plots
xrd_geometry_plot <- xrd_tau_plot +
  xrd_plane_angle

# Front and back strain
# Back strain S4 calculation
xrd_s4_df <- xrd_angles %>%
  dplyr::select(-c(CPdP.CPdN_plane_angle,
                   C.Pd.N_angle,
                   P_transC.Pd.N_angle,
                   P_transC.Pd.C_angle,
                   P_transN.Pd.P_transC_angle,
                   P_transN.Pd.C_angle)) %>%
  mutate(S4_prime_P_trans_N = (Pd.P_transN.C1_angle
    + Pd.P_transN.C2_angle
    + Pd.P_transN.C3.bridge._angle)
    -
    (C1.P_transN.C2_angle
    + C1.P_transN.C3.bridge._angle
    + C2.P_transN.C3.bridge._angle),
    S4_prime_P_trans_C = (Pd.P_transC.C3.bridge._angle
           + Pd.P_transC.C2_angle
           + Pd.P_transC.C1_angle) 
         - (C2.P_transC.C3.bridge._angle 
            + C1.P_transC.C3.bridge._angle
            + C1.P_transC.C2_angle))

# S4' plot for phosphorous trans to N
xrd_s4_plot_ptN <- ggplot(data = xrd_s4_df) +
  geom_label_repel(data = . %>%
                     filter(Structure_source == "New"),
                   aes(x = Structure_source,
                       y = S4_prime_P_trans_N,
                       label = Code),
                   nudge_y = 1) +
  geom_point(aes(x = Structure_source,
                 y = S4_prime_P_trans_N,
                 color = as.character(Coordination_Number)),
             size = 4) +
  scale_color_manual(values = c("#636C9D",
                                "#B71B1B")) +
  labs(y = expression(paste({P}^{trans-N},
                " Symmetric Deformation Coordinate (S4')")),
       x = "Structure") +
  themething +
  theme(legend.position = "none")

# S4' plot for phosphorous trans to C
xrd_s4_plot_ptC <- ggplot(data = xrd_s4_df) +
  geom_label_repel(data = . %>%
                     filter(Structure_source == "New"),
                   aes(x = Structure_source,
                       y = S4_prime_P_trans_C,
                       label = Code),
                   nudge_y = 1) +
  geom_point(aes(x = Structure_source,
                 y = S4_prime_P_trans_C,
                 color = as.character(Coordination_Number)),
             size = 4) +
  scale_color_manual(values = c("#636C9D",
                                "#B71B1B")) +
  labs(y = expression(paste({P}^{trans-C},
                            " Symmetric Deformation Coordinate (S4')")),
       x = "Structure") +
  themething +
  theme(legend.position = "none")

# Combine S4 plots
xrd_s4_plot <- xrd_s4_plot_ptN +
  xrd_s4_plot_ptC

# Front strain CPPd angle differences
xrd_front_strain_df <- xrd_s4_df %>%
  mutate(theo_angle_PtN_Pd = Pd.P_transN.C1_angle
                             + Pd.P_transN.C2_angle
                             + Pd.P_transN.C3.bridge._angle,
         theo_angle_PtC_Pd = Pd.P_transC.C1_angle
         + Pd.P_transC.C2_angle
         + Pd.P_transC.C3.bridge._angle) %>%
  rowwise() %>%
  mutate(max_dev_angle_PtN = sort(c((Pd.P_transN.C1_angle 
                             - theo_angle_PtN_Pd/3),
                             (Pd.P_transN.C2_angle 
                             - theo_angle_PtN_Pd/3),
                             (Pd.P_transN.C3.bridge._angle 
                             - theo_angle_PtN_Pd/3)),
                             decreasing = T)[1],
         max_dev_angle_PtC = sort(c((Pd.P_transC.C1_angle 
                                     - theo_angle_PtC_Pd/3),
                                    (Pd.P_transC.C2_angle 
                                     - theo_angle_PtC_Pd/3),
                                    (Pd.P_transC.C3.bridge._angle 
                                     - theo_angle_PtC_Pd/3)),
                                  decreasing = T)[1]) %>%
  ungroup()

# Trisector deviation plot P trans N
xrd_front_strain_plot_PtN <- ggplot(data = xrd_front_strain_df) +
  geom_label_repel(data = . %>%
                     filter(Structure_source == "New"),
                   aes(x = Structure_source,
                       y = max_dev_angle_PtN,
                       label = Code),
                   nudge_y = 1) +
  geom_point(aes(x = Structure_source,
                 y = max_dev_angle_PtN,
                 color = as.character(Coordination_Number)),
             size = 4) +
  scale_color_manual(values = c("#636C9D",
                                "#B71B1B")) +
  labs(y = expression(paste("Largest deviation from no-tilt ",
                                P^{trans-N},
                                "-Pd bond (\u00B0)")),
       x = "Structure") +
  themething +
  theme(legend.position = "none")

# Trisector deviation plot P trans N
xrd_front_strain_plot_PtC <- ggplot(data = xrd_front_strain_df) +
  geom_label_repel(data = . %>%
                     filter(Structure_source == "New"),
                   aes(x = Structure_source,
                       y = max_dev_angle_PtC,
                       label = Code),
                   nudge_y = 1) +
  geom_point(aes(x = Structure_source,
                 y = max_dev_angle_PtC,
                 color = as.character(Coordination_Number)),
             size = 4) +
  scale_color_manual(values = c("#636C9D",
                                "#B71B1B")) +
  labs(y = expression(paste("Largest Deviation from no-tilt ",
                            P^{trans-C},
                            "-Pd bond (\u00B0)")),
       x = "Structure") +
  themething +
  theme(legend.position = "none")

# Combine front strain plots
xrd_front_strain_plot <- xrd_front_strain_plot_PtN +
  xrd_front_strain_plot_PtC

#===========
# Manuscript
#===========
# Construct df
fs_df_ms <- xrd_s4_df %>%
  left_join(xrd_front_strain_df %>%
              dplyr::select(
                max_dev_angle_PtN,
                Name),
            by = "Name",
            relationship = "many-to-many") %>%
  filter(Structure_source == "New"
         & Name != "F[4]*-N^{Cl}*Ms-OMe-P^t*Bu[3]-MeCN")

# Plot for Fig. 2
back_front_strain_plot <- ggplot() +
  geom_ribbon(data = data.frame("Potential" = c(-5, 
                                                5),
                                "Upper" = c(45, 
                                            45),
                                "Lower" = c(28, 
                                            28)),
              aes(x = Potential,
                  ymax = Upper,
                  ymin = Lower),
              fill = "#B71B1B",
              alpha = 0.2) +
  annotate(geom = "text",
           x = 1,
           y = 25,
           label = "Low front and \nback strain",
           size = 3.5,
           color = "#B71B1B") +
  geom_ribbon(data = data.frame("Potential" = c(10, 
                                                30),
                                "Upper" = c(10, 
                                            10),
                                "Lower" = c(-15, 
                                            -15)),
              aes(x = Potential,
                  ymax = Upper,
                  ymin = Lower),
              fill = "#636C9D",
              alpha = 0.2) +
  annotate(geom = "text",
           x = 24,
           y = 13.3,
           label = "High front and \nback strain",
           size = 3.5,
           color = "#636C9D") +
  geom_label_repel(data = fs_df_ms %>%
                     filter(S4_prime_P_trans_N > 24),
                   aes(
                     x = max_dev_angle_PtN,
                     y = S4_prime_P_trans_N,
                     label = Name),
                   parse = T,
                   nudge_y = 8,
                   nudge_x = 5,
                   size = 3.5) +
  geom_label_repel(data = fs_df_ms %>%
                     filter(S4_prime_P_trans_N < 24
                            & S4_prime_P_trans_N > 0),
                   aes(
                     x = max_dev_angle_PtN,
                     y = S4_prime_P_trans_N,
                     label = Name),
                   parse = T,
                   nudge_y = -5,
                   nudge_x = 5,
                   size = 3.5) +
  geom_label_repel(data = fs_df_ms %>%
                     filter(Coordination_Number == 3),
                   aes(
                     x = max_dev_angle_PtN,
                     y = S4_prime_P_trans_N,
                     label = Name),
                   parse = T,
                   nudge_y = 7,
                   nudge_x = 3,
                   size = 3.5) + 
  geom_point(data = fs_df_ms,
             aes(fill = 
                   as.factor(
                     Coordination_Number),
                 x = max_dev_angle_PtN,
                 y = S4_prime_P_trans_N),
             pch = 21,
             size = 4) +
  scale_fill_manual(values = c("#636C9D",
                               "#B71B1B")) +
  coord_cartesian(xlim = c(0,25),
                  ylim = c(-10,40)) +
  labs(fill = "Coord. No.",
       x = "Trisector deviation (φ, °)",
       y = "S4' P(trans-N)") +
  theme(legend.position = "none") +
  themething
