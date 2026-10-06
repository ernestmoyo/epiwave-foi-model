# ==============================================================================
# Generate animated visualisations for the EpiWave FOI presentation
# Output: GIFs in presentations/slidev/public/images/
# This script is SEPARATE from the main model — presentation only.
# ==============================================================================

library(deSolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gganimate)
library(gifski)

project_root <- "c:/Users/ernes/Documents/R/epiwave-foi-model"
img_dir <- file.path(project_root, "presentations/slidev/public/images")
dir.create(img_dir, showWarnings = FALSE, recursive = TRUE)

# Source the model for Stage 1 functions only (no greta needed)
source(file.path(project_root, "R/epiwave-foi-model.R"))

# ==============================================================================
# Animation 1: ODE dynamics building up over time
# Shows how x (human prev) and z (mosquito prev) evolve at one site
# ==============================================================================
message("Generating: ODE dynamics animation...")

n_sites <- 1; n_times <- 48
times <- seq(0, n_times * 30, by = 30)

m <- get_fixed_m(times, "Site_01", baseline_m = 2.0, seasonal_amplitude = 0.6)
a <- get_fixed_a(times, "Site_01", baseline_a = 0.3)
g <- get_fixed_g(times, "Site_01", baseline_g = 1/10)

itn <- matrix(seq(0, 0.7, length.out = length(times)), ncol = 1)
adj <- apply_interventions(m, a, g, itn_coverage = itn, itn_susceptibility = 0.8)

ode_sol <- solve_ross_macdonald_multi_site(adj$m, adj$a, adj$g, times)

ode_df <- data.frame(
  month = rep(seq_along(times), 2),
  value = c(ode_sol$x[, 1], ode_sol$z[, 1]),
  state = rep(c("x (human prevalence)", "z (mosquito prevalence)"),
              each = length(times))
)

p_ode <- ggplot(ode_df, aes(x = month, y = value, colour = state)) +
  geom_line(linewidth = 1.2) +
  scale_colour_manual(values = c("#2E75B6", "#C00000")) +
  labs(title = "Ross-Macdonald ODE Dynamics",
       subtitle = "Month: {frame_along}",
       x = "Month", y = "Prevalence", colour = "") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold")) +
  transition_reveal(month)

anim_save(file.path(img_dir, "anim_ode_dynamics.gif"),
          animate(p_ode, nframes = 100, fps = 15,
                  width = 700, height = 400, renderer = gifski_renderer()))
message("  Saved: anim_ode_dynamics.gif")


# ==============================================================================
# Animation 2: Mechanistic I* vs observed cases building up
# Shows the two-stage concept: I* provides the shape, cases are noisy samples
# ==============================================================================
message("Generating: I* vs cases animation...")

I_star <- compute_mechanistic_prediction(adj$m, adj$a, 0.8, ode_sol$z)

set.seed(456)
reporting_rate <- 0.1
pop <- 10000
true_incidence <- I_star[, 1] * pop
observed_cases <- rpois(length(times), reporting_rate * true_incidence)

stage1_df <- data.frame(
  month = seq_along(times),
  mechanistic = true_incidence,
  cases = observed_cases
)

p_stage1 <- ggplot(stage1_df, aes(x = month)) +
  geom_line(aes(y = mechanistic, colour = "Mechanistic I* x N"),
            linewidth = 1.2, linetype = "dashed") +
  geom_point(aes(y = cases, colour = "Observed Cases"),
             size = 2.5, alpha = 0.7) +
  scale_colour_manual(values = c("Mechanistic I* x N" = "#C00000",
                                 "Observed Cases" = "black")) +
  labs(title = "Stage 1: Mechanistic Prediction vs Observed",
       subtitle = "Month: {frame_along}  |  Reporting rate = 10%",
       x = "Month", y = "Incidence / Cases", colour = "") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold")) +
  transition_reveal(month)

anim_save(file.path(img_dir, "anim_istar_vs_cases.gif"),
          animate(p_stage1, nframes = 100, fps = 15,
                  width = 700, height = 400, renderer = gifski_renderer()))
message("  Saved: anim_istar_vs_cases.gif")


# ==============================================================================
# Animation 3: GP residuals spreading across sites over time
# Heatmap that fills in month by month
# ==============================================================================
message("Generating: GP residuals animation...")

n_sites <- 10
set.seed(123)
sp_coords <- matrix(runif(n_sites * 2, -5, 5), ncol = 2)
sp_norm <- cbind(
  (sp_coords[, 1] - min(sp_coords[, 1])) / diff(range(sp_coords[, 1])),
  (sp_coords[, 2] - min(sp_coords[, 2])) / diff(range(sp_coords[, 2]))
)

set.seed(321)
eps_mat <- simulate_gp_residuals(sp_norm, n_times = 24, sigma = 0.6,
                                  phi = 3.0, rho = 0.75)

eps_df <- expand.grid(site = seq_len(n_sites), month = 1:24)
eps_df$epsilon <- as.vector(eps_mat)

p_gp <- ggplot(eps_df, aes(x = month, y = factor(site), fill = epsilon)) +
  geom_tile() +
  scale_fill_gradient2(low = "#C00000", mid = "white", high = "#2E75B6",
                       midpoint = 0, limits = c(-2, 2)) +
  labs(title = "GP Residuals: Spatial Correlation + AR(1) Temporal",
       subtitle = "Month: {frame_time}",
       x = "Month", y = "Site", fill = "epsilon") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold")) +
  transition_time(month) +
  shadow_mark(past = TRUE)

anim_save(file.path(img_dir, "anim_gp_residuals.gif"),
          animate(p_gp, nframes = 48, fps = 4,
                  width = 700, height = 400, renderer = gifski_renderer()))
message("  Saved: anim_gp_residuals.gif")


# ==============================================================================
# Animation 4: WITH offset vs WITHOUT (I*=0) prediction building up
# The key comparison Nick specified
# ==============================================================================
message("Generating: WITH vs WITHOUT comparison animation...")

n_times_comp <- 48
times_comp <- seq(0, n_times_comp * 30, by = 30)

m_c <- get_fixed_m(times_comp, "Site_01", baseline_m = 2.0, seasonal_amplitude = 0.6)
a_c <- get_fixed_a(times_comp, "Site_01", baseline_a = 0.3)
g_c <- get_fixed_g(times_comp, "Site_01", baseline_g = 1/10)
itn_c <- matrix(seq(0, 0.7, length.out = length(times_comp)), ncol = 1)
adj_c <- apply_interventions(m_c, a_c, g_c, itn_coverage = itn_c)
ode_c <- solve_ross_macdonald_multi_site(adj_c$m, adj_c$a, adj_c$g, times_comp)
I_star_c <- compute_mechanistic_prediction(adj_c$m, adj_c$a, 0.8, ode_c$z)

set.seed(321)
sp1 <- matrix(c(0.5, 0.5), ncol = 2)
eps_1site <- simulate_gp_residuals(sp1, length(times_comp), 0.6, 3.0, 0.75)

pop_c <- 10000
I_true <- exp(0 + log(pmax(I_star_c[, 1], 1e-6)) + as.vector(eps_1site))
set.seed(456)
obs_c <- rpois(length(times_comp), 0.1 * I_true * pop_c)

# WITH: uses I* shape
pred_with <- 0.1 * I_star_c[, 1] * pop_c
# WITHOUT: flat (grand mean)
pred_without <- rep(mean(obs_c), length(times_comp))

comp_df <- data.frame(
  month = rep(seq_along(times_comp), 3),
  value = c(obs_c, pred_with, pred_without),
  model = rep(c("Observed Cases", "WITH I* offset", "WITHOUT (I*=0)"),
              each = length(times_comp))
)

p_comp <- ggplot(comp_df, aes(x = month, y = value, colour = model)) +
  geom_point(data = filter(comp_df, model == "Observed Cases"),
             size = 2, alpha = 0.5) +
  geom_line(data = filter(comp_df, model != "Observed Cases"),
            linewidth = 1.2) +
  scale_colour_manual(values = c("Observed Cases" = "black",
                                 "WITH I* offset" = "#2E75B6",
                                 "WITHOUT (I*=0)" = "#C00000")) +
  labs(title = "WITH Mechanistic Offset vs WITHOUT (I*=0)",
       subtitle = "Month: {frame_along}",
       x = "Month", y = "Expected Cases", colour = "") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold")) +
  transition_reveal(month)

anim_save(file.path(img_dir, "anim_with_vs_without.gif"),
          animate(p_comp, nframes = 100, fps = 15,
                  width = 700, height = 400, renderer = gifski_renderer()))
message("  Saved: anim_with_vs_without.gif")


# ==============================================================================
# Animation 5: Full pipeline flow — data moving through stages
# ==============================================================================
message("Generating: Pipeline flow animation...")

pipeline_df <- data.frame(
  step = 1:6,
  stage = c("Vector Atlas\nParameters",
            "Ross-Macdonald\nODE",
            "Mechanistic\nI* = m*a*b*z",
            "GP Residuals\nε ~ GP(0, K)",
            "Dual Likelihood\nPoisson + Binomial",
            "Posterior\nPredictions"),
  x = c(1, 2, 3, 4, 5, 6),
  y = c(0, 0, 0, 0, 0, 0),
  colour = c("Stage 1", "Stage 1", "Stage 1",
             "Stage 2", "Stage 2", "Stage 2")
)

p_pipe <- ggplot(pipeline_df, aes(x = x, y = y)) +
  geom_point(aes(colour = colour), size = 20, alpha = 0.8) +
  geom_text(aes(label = stage), size = 3, fontface = "bold") +
  geom_segment(data = data.frame(x = 1:5, xend = 2:6, y = 0, yend = 0),
               aes(x = x + 0.3, xend = xend - 0.3, y = y, yend = yend),
               arrow = arrow(length = unit(0.2, "cm")),
               colour = "grey40", linewidth = 1) +
  scale_colour_manual(values = c("Stage 1" = "#2E75B6", "Stage 2" = "#27AE60")) +
  xlim(0.3, 6.7) + ylim(-0.5, 0.5) +
  labs(title = "EpiWave FOI Model Pipeline") +
  theme_void(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5)) +
  transition_manual(step, cumulative = TRUE) +
  enter_fade()

anim_save(file.path(img_dir, "anim_pipeline.gif"),
          animate(p_pipe, nframes = 12, fps = 2,
                  width = 800, height = 250, renderer = gifski_renderer()))
message("  Saved: anim_pipeline.gif")

message("\n=== All animations generated in: ", img_dir, " ===")
