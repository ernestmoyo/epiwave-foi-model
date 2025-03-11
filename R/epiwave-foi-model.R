############################################################
# EpiWave_FOI_Model.R
# A demonstration script illustrating the EpiWave FOI workflow
# with a simple greta model for Bayesian inference
# Author: Ernest Moyo
# Date: 2025-03-11
############################################################

############################################################
# 0. Load Required Libraries
############################################################
library(deSolve)    # for ODE solutions
library(dplyr)      # for data wrangling
library(sp)         # for spatial objects
library(gstat)      # for variogram/kriging
library(ggplot2)    # for basic plotting
library(greta)      # for Bayesian modeling


############################################################
# 1. Input Data: Synthetic Example
############################################################

# 1a. Set seed for reproducibility
set.seed(123)

# 1b. Create a small grid of locations (latitude & longitude)
n_locations <- 30
locs <- data.frame(
  id  = 1:n_locations,
  lon = runif(n_locations, 30, 35),  # i.e., East Africa range since I will be using Tanzania in my work
  lat = runif(n_locations, -5,  0)
)

# 1c. Synthetic spatial covariates:
#     - population
#     - vector_index (proxy for mosquito abundance)
#     - ITN_coverage (proxy for intervention coverage)
locs$population   <- rpois(n_locations, lambda = 1000)
locs$vector_index <- runif(n_locations, 0.1, 1.0)
locs$ITN_coverage <- runif(n_locations, 0, 1)

# 1d. Synthetic disease data:
#     - observed_cases
#     - observed_prevalence
locs$observed_cases      <- rpois(n_locations, locs$population * 0.01)
locs$observed_prevalence <- runif(n_locations, 0.00, 0.20)

# Show the first few rows
cat("Synthetic data (first few rows):\n")
print(head(locs))


############################################################
# 2. Ross-Macdonald Mechanistic Model
############################################################

# The Ross-Macdonald model describes malaria transmission
# between humans (x) and mosquitoes (z):
#
# dx/dt = m * a * b * z * (1 - x) - r * x
# dz/dt = a * c * x * (1 - z) - g * z
#
# We'll demonstrate solving it with deSolve for a single location
# with fixed parameters.

ross_macdonald <- function(t, state, parameters) {
  x <- state[1]  # human infection prevalence
  z <- state[2]  # mosquito infection prevalence

  with(as.list(parameters), {
    dx <- m*a*b*z*(1 - x) - r*x
    dz <- a*c*x*(1 - z) - g*z
    list(c(dx, dz))
  })
}

# Example parameters
params <- c(
  m = 10,   # mosquito density per human
  a = 0.3,  # biting rate
  b = 0.5,  # prob(mosquito -> human infection)
  c = 0.5,  # prob(human -> mosquito infection)
  r = 0.1,  # human recovery rate
  g = 0.2   # mosquito mortality rate
)

# Initial conditions
init_state <- c(x = 0.01, z = 0.01)

# Time range
times <- seq(0, 100, by = 1)

# Solve ODE
out_ODE <- ode(y = init_state, times = times, func = ross_macdonald, parms = params)
out_ODE_df <- as.data.frame(out_ODE)

cat("\nRoss-Macdonald ODE (single example) solved over time.\n")
print(head(out_ODE_df, 10))

# Quick plot
ggplot(out_ODE_df, aes(x = time)) +
  geom_line(aes(y = x, color = "Human (x)")) +
  geom_line(aes(y = z, color = "Mosquito (z)")) +
  labs(x = "Time", y = "Prevalence", title = "Ross-Macdonald ODE Example") +
  theme_minimal() +
  guides(color = guide_legend(title = "State"))


############################################################
# 3. Map Spatial Covariates to Transmission Parameters
############################################################

# NOTE TO SELCH # CHECK THIS LATER
####Note this is to be done when working on a full model##
# When working on a full model, I will calibrate or model how
# covariates (vector_index, ITN_coverage, etc.) inform
# the Ross-Macdonald parameters (m, a, r, g, etc.).
# Here, we use simple hypothetical formulas:

map_parameters <- function(vector_index, itn) {
  m_est <- 5 + 20 * vector_index  # grows with vector index
  a_est <- 0.4 * (1 - itn)        # biting rate decreases with ITN
  r_est <- 0.1                    # kept constant
  g_est <- 0.2 + 0.3 * itn        # higher ITN -> higher mosquito mortality
  list(m = m_est, a = a_est, b = 0.5, c = 0.5, r = r_est, g = g_est)
}

locs <- locs %>%
  rowwise() %>%
  mutate(
    param_list = list(map_parameters(vector_index, ITN_coverage)),
    m = param_list$m,
    a = param_list$a,
    b = param_list$b,
    c = param_list$c,
    r = param_list$r,
    g = param_list$g
  ) %>%
  ungroup()

cat("\nParameters mapped from covariates (first few rows):\n")
print(head(locs))


############################################################
# 4. Deriving FOI (Force of Infection)
############################################################

# FOI_{l,t} = m * a * b * z_{l,t}
# We need z_{l,t}, the mosquito infection prevalence.
# We'll assume an equilibrium solution for demonstration.

library(rootSolve)

equilibrium_prez <- function(params) {
  # We solve the steady-state:
  # 0 = m*a*b*z*(1 - x) - r*x
  # 0 = a*c*x*(1 - z) - g*z

  eq_fn <- function(vars) {
    x <- vars[1]
    z <- vars[2]
    with(as.list(params), {
      eq1 <- m*a*b*z*(1 - x) - r*x
      eq2 <- a*c*x*(1 - z) - g*z
      c(eq1, eq2)
    })
  }

  # Use multiroot with random starts
  start_vals <- c(runif(1, 0, 0.1), runif(1, 0, 0.1))
  sol <- multiroot(eq_fn, start = start_vals)
  sol$root
}

# Compute equilibrium for each location
locs_equil <- locs %>%
  rowwise() %>%
  mutate(
    eq = list(equilibrium_prez(
      params = c(m = m, a = a, b = b, c = c, r = r, g = g)
    )),
    x_eq = eq[1],
    z_eq = eq[2],
    FOI = m * a * b * eq[2]  # FOI = m*a*b*z_eq
  ) %>%
  ungroup()

cat("\nEquilibrium solutions and FOI (first few rows):\n")
print(head(locs_equil))


############################################################
# 5. Refine FOI with a Gaussian Process
############################################################

# We'll treat log(FOI) as a spatial field and use kriging
# to produce a smoothed estimate.

coordinates(locs_equil) <- ~lon+lat
locs_equil$log_FOI <- log(pmax(locs_equil$FOI, 1e-9))

# Empirical variogram
vg <- variogram(log_FOI ~ 1, data = locs_equil)
vg_fit <- fit.variogram(vg, model = vgm("Exp"))

# Ordinary kriging (here, just at the same locations for demo)
krig_res <- krige(
  formula   = log_FOI ~ 1,
  locations = locs_equil,
  newdata   = locs_equil,
  model     = vg_fit
)

# Merge the kriged results
locs_equil$log_FOI_refined <- krig_res$var1.pred
locs_equil$FOI_refined     <- exp(krig_res$var1.pred)

# Convert back to data frame for further steps
locs_equil_df <- as.data.frame(locs_equil)

cat("\nRefined FOI via GP (kriging):\n")
print(head(locs_equil_df))


############################################################
# 6. Clinical Observations & Prevalence
############################################################

# Suppose we want to compare predicted cases to observed cases.
# We'll do a crude predicted_cases measure:
locs_equil_df <- locs_equil_df %>%
  mutate(
    predicted_cases = population * (x_eq * 0.1)  # toy factor: 10% of infected are clinical
  )

comparison_df <- data.frame(
  id        = locs_equil_df$id,
  observed  = locs_equil_df$observed_cases,
  predicted = locs_equil_df$predicted_cases
)

cat("\nObserved vs. predicted clinical cases (example model):\n")
print(head(comparison_df))


############################################################
# 7. Bayesian Inference with greta
############################################################

# NOTE TO SELCH # CHECK THIS LATER
####Note this is to be done when working on a full model##
# In a real FOI model, will incorporate the Ross-Macdonald
# equations (or equilibrium) and the GP into a single
# hierarchical model.
#
# Below is a minimal example of using greta to estimate a
# parameter from the observed case data. We illustrate how
# to set up greta distributions and sample from the posterior.
#
# We'll do a example model:
#   observed_cases[i] ~ Poisson(lambda * population[i])
#   lambda ~ half-normal(0, 10)
#
# Can replace this with more sophisticated models
# (e.g. linking lambda to FOI or x_eq, etc.).

# Prepare observed data
pop <- locs_equil_df$population
obs_cases <- locs_equil_df$observed_cases

# Define greta parameters (with priors)
lambda <- normal(0, 10, truncation = c(0, Inf))  # half-normal

# Mean of Poisson depends on population and lambda
mean_cases <- lambda * pop

# Define likelihood
distribution(obs_cases) <- poisson(mean_cases)

# Define the greta model
m <- model(lambda)

# Sample from the posterior
draws <- mcmc(m, n_samples = 2000, warmup = 500, chains = 2)

cat("\nSummary of greta MCMC for 'lambda':\n")
print(summary(draws))

# Print a quick traceplot for lambda
plot(draws)


############################################################
# 8. Model Outputs & Predictive Maps
############################################################

# Could now use posterior draws for lambda (or other
# parameters) in combination with your FOI or Ross-Macdonald
# equations to generate full predictive maps. Below is an
# example map of the refined FOI in the location grid.

ggplot(locs_equil_df, aes(x = lon, y = lat)) +
  geom_point(aes(size = FOI_refined)) +
  labs(
    title = "Refined FOI (Gaussian Process)",
    x = "Longitude",
    y = "Latitude",
    size = "FOI"
  ) +
  theme_minimal()

