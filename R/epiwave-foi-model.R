# ==============================================================================
# EpiWave FOI Model: Vector-Informed Malaria Transmission Mapping
# ==============================================================================
# Author: Ernest Moyo (NM-AIST / Vector Atlas)
# PhD Objective 2 — Framework Development
#
# Two-stage computational strategy:
#   Stage 1: Fixed entomological params (m,a,g) -> Ross-Macdonald ODE -> I*(s,t)
#   Stage 2: GP residuals + dual likelihood (Poisson cases + Binomial prevalence)
#            with log(I*) as fixed offset
#
# Key design choices (per supervisor Prof Nick Golding):
#   - Do NOT infer dynamic ODE parameters. Use fixed Vector Atlas / temperature
#     model estimates, solve ODEs ONCE per pixel, then calibrate to case data.
#   - Use a GP for spatially-correlated residuals (NOT NegBin overdispersion).
#   - Dual likelihood (cases + prevalence) makes alpha and gamma identifiable.
#
# Reference: https://github.com/idem-lab/epiwave.mapping
# ==============================================================================

library(deSolve)
library(dplyr)
library(tidyr)
library(ggplot2)
# greta + greta.gp loaded lazily inside fit_epiwave_gp()

# ==============================================================================
# STAGE 1: FIXED ENTOMOLOGICAL PARAMETERS
# ==============================================================================

#' Generate mosquito-to-human ratio (m) with bimodal East African seasonality.
#' @param time Numeric vector of time points (days)
#' @param location Character vector of location IDs
#' @param vector_atlas_data Optional data.frame (time, location, m) for empirical values
#' @param baseline_m Baseline ratio (default 2.0)
#' @param seasonal_amplitude Seasonal amplitude [0,1] (default 0.5)
#' @param phase_shift Phase shift in radians (default 0)
#' @return Matrix [n_times x n_locations]
#' @export
get_fixed_m <- function(time, location, vector_atlas_data = NULL,
                        baseline_m = 2.0, seasonal_amplitude = 0.5,
                        phase_shift = 0) {
  n_times <- length(time)
  n_locs  <- length(location)

  if (!is.null(vector_atlas_data)) {
    m_matrix <- vector_atlas_data %>%
      select(time, location, m) %>%
      pivot_wider(names_from = location, values_from = m) %>%
      select(-time) %>%
      as.matrix()
    if (!all(dim(m_matrix) == c(n_times, n_locs)))
      stop("Vector Atlas data dimensions do not match time/location specifications")
    return(m_matrix)
  }

  m_matrix <- matrix(NA, nrow = n_times, ncol = n_locs)
  for (i in 1:n_locs) {
    year <- time / 365
    seasonal_factor <- 1 + seasonal_amplitude * (
      sin(2 * 2 * pi * year + phase_shift) - 0.3 * cos(2 * pi * year)
    )
    m_matrix[, i] <- baseline_m * pmax(seasonal_factor, 0.1)
  }
  return(m_matrix)
}


#' Generate human biting rate (a), optionally temperature-dependent.
#' @param time Numeric vector of time points
#' @param location Character vector of location IDs
#' @param temperature_data Optional data.frame or matrix of temperatures
#' @param baseline_a Baseline biting rate per mosquito per day (default 0.3)
#' @return Matrix [n_times x n_locations]
#' @export
get_fixed_a <- function(time, location, temperature_data = NULL,
                        baseline_a = 0.3) {
  n_times <- length(time)
  n_locs  <- length(location)

  if (is.null(temperature_data)) {
    return(matrix(baseline_a, nrow = n_times, ncol = n_locs))
  }

  a_matrix <- matrix(NA, nrow = n_times, ncol = n_locs)
  if (is.data.frame(temperature_data)) {
    for (i in 1:n_locs) {
      loc_data <- temperature_data[temperature_data$location == location[i], ]
      if (nrow(loc_data) > 0) {
        temp_effect <- pmax(1 - 0.02 * (loc_data$temperature - 25)^2, 0.2)
        a_matrix[, i] <- baseline_a * temp_effect
      } else {
        a_matrix[, i] <- baseline_a
      }
    }
  } else {
    a_matrix <- baseline_a * pmax(1 - 0.02 * (temperature_data - 25)^2, 0.2)
  }
  return(a_matrix)
}


#' Generate mosquito mortality rate (g), optionally temperature-dependent.
#' @param time Numeric vector of time points
#' @param location Character vector of location IDs
#' @param temperature_data Optional data.frame or matrix of temperatures
#' @param humidity_data Optional (reserved for future use)
#' @param baseline_g Baseline mortality rate per day (default 1/10)
#' @return Matrix [n_times x n_locations]
#' @export
get_fixed_g <- function(time, location, temperature_data = NULL,
                        humidity_data = NULL, baseline_g = 1/10) {
  n_times <- length(time)
  n_locs  <- length(location)

  if (is.null(temperature_data)) {
    return(matrix(baseline_g, nrow = n_times, ncol = n_locs))
  }

  g_matrix <- matrix(NA, nrow = n_times, ncol = n_locs)
  if (is.data.frame(temperature_data)) {
    for (i in 1:n_locs) {
      loc_data <- temperature_data[temperature_data$location == location[i], ]
      if (nrow(loc_data) > 0) {
        temp_effect <- 1 + 0.01 * (loc_data$temperature - 25)^2
        g_matrix[, i] <- baseline_g * temp_effect
      } else {
        g_matrix[, i] <- baseline_g
      }
    }
  } else {
    g_matrix <- baseline_g * (1 + 0.01 * (temperature_data - 25)^2)
  }
  return(g_matrix)
}


#' Apply ITN/IRS intervention effects to entomological parameters.
#' Based on Griffin et al. (2010) and Bhatt et al. (2015).
#' @param m,a,g Matrices [n_times x n_locations] of baseline parameters
#' @param itn_coverage Matrix of ITN coverage [0,1] (optional)
#' @param irs_coverage Matrix of IRS coverage [0,1] (optional)
#' @param resistance_index Insecticide resistance [0=susceptible, 1=resistant] (default 0)
#' @param itn_kill_rate Kill rate among ITN encounters (default 0.5)
#' @param itn_feeding_inhibit Feeding reduction from ITNs (default 0.3)
#' @param itn_mortality_boost Mortality increase from ITNs (default 0.3)
#' @param irs_efficacy IRS killing efficacy (default 0.5)
#' @return List with adjusted m, a, g matrices
#' @export
apply_interventions <- function(m, a, g, itn_coverage = NULL, irs_coverage = NULL,
                                resistance_index = 0, itn_kill_rate = 0.5,
                                itn_feeding_inhibit = 0.3, itn_mortality_boost = 0.3,
                                irs_efficacy = 0.5) {
  u <- 1 - resistance_index
  m_adj <- m; a_adj <- a; g_adj <- g

  if (!is.null(itn_coverage)) {
    if (any(itn_coverage < 0 | itn_coverage > 1, na.rm = TRUE))
      stop("ITN coverage values must be in range [0, 1]")
    m_adj <- m_adj * ((1 - itn_coverage) + itn_coverage * (1 - u * itn_kill_rate))
    a_adj <- a_adj * (1 - itn_coverage * u * itn_feeding_inhibit)
    g_adj <- g_adj * (1 + itn_coverage * u * itn_mortality_boost)
  }

  if (!is.null(irs_coverage)) {
    if (any(irs_coverage < 0 | irs_coverage > 1, na.rm = TRUE))
      stop("IRS coverage values must be in range [0, 1]")
    g_adj <- g_adj * (1 + irs_coverage * u * irs_efficacy)
    a_adj <- a_adj * (1 - irs_coverage * u * 0.1)
  }

  list(m = m_adj, a = a_adj, g = g_adj)
}


# ==============================================================================
# STAGE 1: ROSS-MACDONALD ODE SOLVER
# ==============================================================================

#' Ross-Macdonald ODE system: dx/dt, dz/dt for human and mosquito prevalence.
#' @keywords internal
ross_macdonald_ode <- function(t, state, parms) {
  x <- state[1]; z <- state[2]
  m <- if (is.function(parms$m)) parms$m(t) else parms$m
  a <- if (is.function(parms$a)) parms$a(t) else parms$a
  g <- if (is.function(parms$g)) parms$g(t) else parms$g
  b <- parms$b; c <- parms$c; r <- parms$r

  dx_dt <- m * a * b * z * (1 - x) - r * x
  dz_dt <- a * c * x * (1 - z) - g * z
  list(c(dx_dt, dz_dt))
}


#' Solve Ross-Macdonald ODEs for multiple sites (one-time forward integration).
#' @param m_matrix,a_matrix,g_matrix Matrices [n_times x n_sites]
#' @param times Numeric vector of time points (days)
#' @param b Transmission prob mosquito->human (default 0.8)
#' @param c Transmission prob human->mosquito (default 0.8)
#' @param r Human recovery rate per day (default 1/7)
#' @param x0,z0 Initial prevalences (default 0.01, 0.001)
#' @return List with x and z matrices [n_times x n_sites]
#' @export
solve_ross_macdonald_multi_site <- function(m_matrix, a_matrix, g_matrix, times,
                                            b = 0.8, c = 0.8, r = 1/7,
                                            x0 = 0.01, z0 = 0.001) {
  n_times <- length(times)
  n_sites <- ncol(m_matrix)

  if (!all(dim(m_matrix) == dim(a_matrix), dim(m_matrix) == dim(g_matrix)))
    stop("Parameter matrices (m, a, g) must have identical dimensions")
  if (nrow(m_matrix) != n_times)
    stop("Number of rows in parameter matrices must match length of times vector")

  x_matrix <- matrix(NA, nrow = n_times, ncol = n_sites)
  z_matrix <- matrix(NA, nrow = n_times, ncol = n_sites)

  for (site in 1:n_sites) {
    m_func <- approxfun(times, m_matrix[, site], rule = 2)
    a_func <- approxfun(times, a_matrix[, site], rule = 2)
    g_func <- approxfun(times, g_matrix[, site], rule = 2)

    parms <- list(m = m_func, a = a_func, g = g_func, b = b, c = c, r = r)
    solution <- ode(y = c(x = x0, z = z0), times = times,
                    func = ross_macdonald_ode, parms = parms, method = "lsoda")
    x_matrix[, site] <- solution[, "x"]
    z_matrix[, site] <- solution[, "z"]
  }

  list(x = x_matrix, z = z_matrix)
}


#' Compute mechanistic incidence prediction: I*(t,s) = m*a*b*z * N
#' @param m_matrix,a_matrix,z_matrix,population_matrix Matrices [n_times x n_sites]
#' @param b Transmission probability mosquito->human
#' @return Matrix I_star [n_times x n_sites]
#' @export
compute_mechanistic_prediction <- function(m_matrix, a_matrix, b,
                                           z_matrix, population_matrix) {
  if (!all(dim(m_matrix) == dim(a_matrix),
           dim(m_matrix) == dim(z_matrix),
           dim(m_matrix) == dim(population_matrix)))
    stop("All input matrices must have identical dimensions")

  m_matrix * a_matrix * b * z_matrix * population_matrix
}


# ==============================================================================
# STAGE 2: BAYESIAN CALIBRATION (GP residuals + dual likelihood)
# ==============================================================================
# Model specification (per Nick Golding / epiwave.mapping):
#   log(I_{l,t}) = alpha + log(I*_{l,t}) + epsilon_{l,t}
#   epsilon ~ GP(0, K), K = sigma^2 * K_space(phi) * K_time(theta)
#   C_{l,t} ~ Poisson(gamma * I_{l,t})       -- case counts
#   Y_{l,t} ~ Binomial(N_{l,t}, x_{l,t})     -- prevalence surveys
# ==============================================================================

#' Build separable space-time GP kernel using greta.gp.
#'
#' K = sigma^2 * Matern52(phi, columns=1:2) * Exponential(theta, columns=3)
#'
#' @param sigma2 GP marginal variance (greta scalar)
#' @param phi Spatial lengthscale (greta scalar)
#' @param theta Temporal lengthscale (greta scalar)
#' @return greta.gp kernel object
#' @export
build_gp_kernel <- function(sigma2, phi, theta) {
  if (!requireNamespace("greta.gp", quietly = TRUE))
    stop("Package 'greta.gp' required. Install with install.packages('greta.gp')")

  # Spatial kernel: isotropic Matern 5/2 as product of 1D kernels per dimension.
  # sigma2 on K_lon only; K_lat has unit variance to avoid double-counting.
  K_lon <- greta.gp::mat52(lengthscales = phi, variance = sigma2, columns = 1L)
  K_lat <- greta.gp::mat52(lengthscales = phi, variance = 1,      columns = 2L)
  K_space <- K_lon * K_lat

  # Temporal kernel: Exponential (Matern 1/2) on column 3 (time)
  K_time <- greta.gp::expo(
    lengthscales = theta,
    variance = 1,
    columns = 3
  )

  # Separable product kernel
  K_space * K_time
}


#' Simulate GP residuals from a separable space-time covariance.
#'
#' Uses MASS::mvrnorm() with the true covariance matrix for data generation.
#'
#' @param coords Matrix [n_obs x 3] of (lon, lat, time_normalised) coordinates
#' @param sigma2 GP marginal variance (scalar)
#' @param phi Spatial lengthscale (scalar)
#' @param theta Temporal lengthscale (scalar)
#' @return Numeric vector of GP residuals (length n_obs)
#' @export
simulate_gp_residuals <- function(coords, sigma2, phi, theta) {
  n <- nrow(coords)

  # Spatial distances (Euclidean on lon, lat)
  spatial_dists <- as.matrix(dist(coords[, 1:2]))
  # Temporal distances
  temporal_dists <- as.matrix(dist(coords[, 3]))

  # Matern 5/2 spatial kernel
  r_s <- sqrt(5) * spatial_dists / phi
  K_space <- sigma2 * (1 + r_s + r_s^2 / 3) * exp(-r_s)

  # Exponential temporal kernel
  K_time <- exp(-temporal_dists / theta)

  # Separable product + jitter for numerical stability
  K <- K_space * K_time + diag(1e-6, n)

  as.vector(MASS::mvrnorm(1, mu = rep(0, n), Sigma = K))
}


#' Simulate prevalence survey data from ODE human prevalence.
#'
#' Generates Binomial survey data at a random subset of (site, time) pairs.
#'
#' @param x_matrix Matrix [n_times x n_sites] of human prevalence from ODE
#' @param gp_adjustment Vector of exp(epsilon) adjustments (length n_times*n_sites)
#' @param survey_fraction Fraction of (site,time) pairs with surveys (default 0.3)
#' @param sample_size_range Range of survey sample sizes (default c(50, 200))
#' @param seed Random seed for reproducibility
#' @return List with survey_indices, n_tested, n_positive, true_prevalence
#' @export
simulate_prevalence_surveys <- function(x_matrix, gp_adjustment = NULL,
                                        survey_fraction = 0.3,
                                        sample_size_range = c(50, 200),
                                        seed = 789) {
  set.seed(seed)
  n_obs <- length(x_matrix)

  # Select random subset of observations to have surveys
  n_surveys <- round(n_obs * survey_fraction)
  survey_indices <- sort(sample(seq_len(n_obs), n_surveys))

  # True prevalence at survey locations (adjusted by GP if provided)
  x_vec <- as.vector(x_matrix)
  true_prev <- x_vec[survey_indices]
  if (!is.null(gp_adjustment)) {
    true_prev <- pmin(true_prev * gp_adjustment[survey_indices], 0.99)
  }
  true_prev <- pmax(true_prev, 1e-6)

  # Survey sample sizes
  n_tested <- sample(sample_size_range[1]:sample_size_range[2],
                     n_surveys, replace = TRUE)

  # Observed positive counts
  n_positive <- rbinom(n_surveys, size = n_tested, prob = true_prev)

  list(
    survey_indices = survey_indices,
    n_tested = n_tested,
    n_positive = n_positive,
    true_prevalence = true_prev
  )
}


#' Fit Stage 2 GP model with dual likelihood via greta.
#'
#' Model:
#'   log(I) = alpha + log(I*) + epsilon, epsilon ~ GP(0, K)
#'   K = sigma^2 * Matern52(phi) * Exponential(theta)
#'   cases ~ Poisson(gamma * I)
#'   prevalence ~ Binomial(N_tested, x_adjusted)
#'
#' @param observed_cases Matrix [n_times x n_sites] of case counts
#' @param I_star Matrix [n_times x n_sites] of Stage 1 mechanistic predictions
#' @param x_star Matrix [n_times x n_sites] of ODE human prevalence
#' @param coords Matrix [n_obs x 3] of (lon, lat, time_normalised) coordinates
#' @param prev_data List from simulate_prevalence_surveys() (NULL for case-only)
#' @param use_mechanistic If TRUE, use log(I*) as offset; if FALSE, intercept-only
#' @param inducing Optional matrix of inducing point coordinates for sparse GP
#' @param gp_tol Jitter for GP numerical stability (default 1e-3)
#' @return greta model object for mcmc()
#' @export
fit_epiwave_gp <- function(observed_cases, I_star, x_star,
                           coords, prev_data = NULL,
                           use_mechanistic = TRUE,
                           inducing = NULL,
                           gp_tol = 1e-3) {
  # --- Input validation ---
  if (!is.matrix(observed_cases))
    stop("observed_cases must be a matrix [n_times x n_sites]")
  if (!identical(dim(observed_cases), dim(I_star)))
    stop("observed_cases and I_star must have identical dimensions")
  if (any(is.na(observed_cases)))
    stop("observed_cases contains NA values")
  if (any(I_star < 0))
    stop("I_star contains negative values")

  if (!requireNamespace("greta", quietly = TRUE))
    stop("Package 'greta' required. Run R/greta_setup.R first.")
  if (!requireNamespace("greta.gp", quietly = TRUE))
    stop("Package 'greta.gp' required.")
  suppressPackageStartupMessages({
    library(greta)
    library(greta.gp)
  })

  I_star_vec <- as.vector(I_star)
  n_obs      <- length(I_star_vec)

  # Replace zeros/negatives in I_star to avoid log(0)
  I_star_clean <- pmax(I_star_vec, 1e-6)
  log_I_star_raw <- log(I_star_clean)
  if (any(!is.finite(log_I_star_raw)))
    stop("log(I_star) contains non-finite values after cleaning")

  # Wrap observed data as greta data arrays
  cases_greta <- as_data(as.vector(observed_cases))

  # --- Priors ---
  alpha     <- normal(0, 1)                                    # intercept
  gamma_rr  <- normal(0.1, 0.05, truncation = c(0.001, Inf))  # reporting rate
  log_sigma <- normal(-1, 1)                                   # log GP variance
  sigma2    <- exp(log_sigma)
  log_phi   <- normal(1, 1)                                    # log spatial lengthscale
  phi       <- exp(log_phi)
  log_theta <- normal(1, 1)                                    # log temporal lengthscale
  theta     <- exp(log_theta)

  # --- GP kernel and residuals ---
  K <- build_gp_kernel(sigma2, phi, theta)

  coord_array <- as_data(coords)
  epsilon <- if (!is.null(inducing)) {
    inducing_array <- as_data(inducing)
    gp(coord_array, K, inducing = inducing_array, tol = gp_tol)
  } else {
    gp(coord_array, K, tol = gp_tol)
  }

  # --- Latent infection incidence ---
  log_I <- if (use_mechanistic) {
    alpha + log_I_star_raw + epsilon
  } else {
    alpha + mean(log_I_star_raw) + epsilon
  }
  I_latent <- exp(log_I)

  # --- Case likelihood (Poisson) ---
  expected_cases <- gamma_rr * I_latent
  distribution(cases_greta) <- poisson(expected_cases)

  # --- Prevalence likelihood (Binomial) if data provided ---
  if (!is.null(prev_data)) {
    x_star_vec <- as.vector(x_star)
    x_at_surveys <- pmax(x_star_vec[prev_data$survey_indices], 1e-6)
    x_at_surveys <- pmin(x_at_surveys, 0.999)

    # GP adjustment to prevalence at survey locations
    # Use ilogit to map adjusted prevalence to valid (0,1) probability
    gp_at_surveys  <- exp(epsilon[prev_data$survey_indices])
    log_odds_base  <- log(x_at_surveys) - log(1 - x_at_surveys)
    log_odds_adj   <- log_odds_base + epsilon[prev_data$survey_indices]
    prev_prob      <- ilogit(log_odds_adj)

    n_pos_greta <- as_data(prev_data$n_positive)
    distribution(n_pos_greta) <- binomial(prev_data$n_tested, prev_prob)
  }

  model(alpha, gamma_rr, log_sigma, log_phi, log_theta)
}


# ==============================================================================
# VALIDATION: SIMULATION-ESTIMATION STUDY
# ==============================================================================

#' Run simulation-estimation study with GP + dual likelihood.
#'
#' Comparison: GP+offset vs GP-only (no mechanistic information).
#'
#' @param n_sites Number of spatial locations (default 10)
#' @param n_times Number of monthly time points (default 48)
#' @param true_params List of true parameter values (NULL for defaults)
#' @param include_interventions Simulate ITN scale-up (default TRUE)
#' @param run_mcmc Whether to run MCMC sampling (default TRUE)
#' @param n_samples,warmup,chains MCMC settings
#' @param use_sparse_gp Use inducing points for sparse GP (default TRUE)
#' @param n_inducing Number of inducing points for sparse GP (default 25)
#' @param draws_gp_with,draws_gp_without Pre-computed draws (optional)
#' @return List of simulation results
#' @export
simulate_and_estimate <- function(n_sites = 10, n_times = 48,
                                  true_params = NULL,
                                  include_interventions = TRUE,
                                  run_mcmc = TRUE,
                                  n_samples = 1000, warmup = 500, chains = 2,
                                  use_sparse_gp = TRUE,
                                  n_inducing = 25,
                                  draws_gp_with = NULL,
                                  draws_gp_without = NULL) {

  # ---- Step 1: True parameters ----
  if (is.null(true_params)) {
    true_params <- list(
      baseline_m = 2.0, baseline_a = 0.3, baseline_g = 1/10,
      b = 0.8, c = 0.8, r = 1/7,
      population = 10000, reporting_rate = 0.1,
      # True GP hyperparameters
      alpha = 0,          # intercept (no systematic bias)
      gp_sigma2 = 0.3,    # GP marginal variance
      gp_phi = 3.0,       # spatial lengthscale
      gp_theta = 0.3      # temporal lengthscale (on normalised time)
    )
  }

  times     <- seq(0, n_times * 30, by = 30)
  locations <- paste0("Site_", sprintf("%02d", 1:n_sites))

  # ---- Step 2: Spatial coordinates ----
  set.seed(123)
  spatial_coords <- matrix(runif(n_sites * 2, min = -5, max = 5),
                           ncol = 2, dimnames = list(locations, c("lon", "lat")))

  # ---- Step 3: Stage 1 — ODE with fixed entomological params ----
  m_true <- get_fixed_m(times, locations, baseline_m = true_params$baseline_m,
                        seasonal_amplitude = 0.6)
  a_true <- get_fixed_a(times, locations, baseline_a = true_params$baseline_a)
  g_true <- get_fixed_g(times, locations, baseline_g = true_params$baseline_g)

  if (include_interventions) {
    itn_coverage <- matrix(rep(seq(0, 0.7, length.out = length(times)), n_sites),
                           nrow = length(times), ncol = n_sites)
    params_adj <- apply_interventions(m = m_true, a = a_true, g = g_true,
                                      itn_coverage = itn_coverage,
                                      resistance_index = 0.2)
    m_true <- params_adj$m; a_true <- params_adj$a; g_true <- params_adj$g
  }

  ode_solution <- solve_ross_macdonald_multi_site(
    m_matrix = m_true, a_matrix = a_true, g_matrix = g_true,
    times = times, b = true_params$b, c = true_params$c, r = true_params$r
  )

  pop_matrix <- matrix(true_params$population, nrow = length(times), ncol = n_sites)

  I_star <- compute_mechanistic_prediction(
    m_matrix = m_true, a_matrix = a_true, b = true_params$b,
    z_matrix = ode_solution$z, population_matrix = pop_matrix
  )

  x_star <- ode_solution$x  # human prevalence from ODE

  # ---- Step 4: Build coordinate matrix (lon, lat, time_normalised) ----
  time_norm <- (times - min(times)) / max(times - min(times))  # [0, 1]
  lon_norm  <- (spatial_coords[, 1] - min(spatial_coords[, 1])) /
    max(spatial_coords[, 1] - min(spatial_coords[, 1]) + 1e-10)
  lat_norm  <- (spatial_coords[, 2] - min(spatial_coords[, 2])) /
    max(spatial_coords[, 2] - min(spatial_coords[, 2]) + 1e-10)

  # Expand to full grid: each row is one (site, time) observation
  coord_grid <- expand.grid(time_idx = seq_along(times), site_idx = seq_len(n_sites))
  coords <- cbind(
    lon  = lon_norm[coord_grid$site_idx],
    lat  = lat_norm[coord_grid$site_idx],
    time = time_norm[coord_grid$time_idx]
  )

  # ---- Step 5: Simulate true GP residuals ----
  set.seed(321)
  epsilon_true <- simulate_gp_residuals(
    coords  = coords,
    sigma2  = true_params$gp_sigma2,
    phi     = true_params$gp_phi,
    theta   = true_params$gp_theta
  )

  # ---- Step 6: Generate observed data with GP residuals ----
  I_star_vec <- as.vector(I_star)
  I_true_vec <- exp(true_params$alpha + log(I_star_vec + 1e-10) + epsilon_true)
  expected_cases <- true_params$reporting_rate * I_true_vec

  set.seed(456)
  observed_cases_vec <- rpois(length(expected_cases), lambda = pmax(expected_cases, 1e-10))
  observed_cases <- matrix(observed_cases_vec, nrow = length(times), ncol = n_sites)

  # Prevalence surveys
  gp_adjustment <- exp(epsilon_true)
  prev_data <- simulate_prevalence_surveys(
    x_matrix = x_star, gp_adjustment = gp_adjustment,
    survey_fraction = 0.3, seed = 789
  )

  message(sprintf("Data generated: %d sites x %d times = %d obs, %d prevalence surveys",
                  n_sites, length(times), length(I_star_vec), length(prev_data$n_positive)))

  # ---- Step 7: Inducing points for sparse GP ----
  inducing <- NULL
  if (use_sparse_gp) {
    n_spatial  <- min(ceiling(sqrt(n_inducing)), n_sites)
    n_temporal <- min(ceiling(n_inducing / n_spatial), length(times))
    ind_s <- seq(1, n_sites, length.out = n_spatial)
    ind_t <- seq(1, length(times), length.out = n_temporal)
    ind_grid <- expand.grid(s = round(ind_s), t = round(ind_t))
    inducing <- cbind(
      lon  = lon_norm[ind_grid$s],
      lat  = lat_norm[ind_grid$s],
      time = time_norm[ind_grid$t]
    )
    message(sprintf("Sparse GP: %d inducing points", nrow(inducing)))
  }

  # ---- Step 8: Fit GP+offset model ----
  if (run_mcmc) {
    message("\n--- Fitting GP + offset model ---")
    model_gp_with <- tryCatch(
      fit_epiwave_gp(
        observed_cases = observed_cases, I_star = I_star, x_star = x_star,
        coords = coords, prev_data = prev_data,
        use_mechanistic = TRUE, inducing = inducing
      ),
      error = function(e) { message("Model build error (GP+offset): ", e$message); NULL }
    )

    if (!is.null(model_gp_with) && is.null(draws_gp_with)) {
      t0 <- proc.time()
      draws_gp_with <- tryCatch(
        greta::mcmc(model_gp_with, n_samples = n_samples,
                    warmup = warmup, chains = chains, verbose = TRUE),
        error = function(e) { message("MCMC error (GP+offset): ", e$message); NULL }
      )
      message(sprintf("GP+offset model: %.0fs", (proc.time() - t0)["elapsed"]))
    }

    # ---- Step 9: Fit GP-only model (no mechanistic offset) ----
    message("\n--- Fitting GP-only model (no offset) ---")
    model_gp_without <- tryCatch(
      fit_epiwave_gp(
        observed_cases = observed_cases, I_star = I_star, x_star = x_star,
        coords = coords, prev_data = prev_data,
        use_mechanistic = FALSE, inducing = inducing
      ),
      error = function(e) { message("Model build error (GP-only): ", e$message); NULL }
    )

    if (!is.null(model_gp_without) && is.null(draws_gp_without)) {
      t0 <- proc.time()
      draws_gp_without <- tryCatch(
        greta::mcmc(model_gp_without, n_samples = n_samples,
                    warmup = warmup, chains = chains, verbose = TRUE),
        error = function(e) { message("MCMC error (GP-only): ", e$message); NULL }
      )
      message(sprintf("GP-only model: %.0fs", (proc.time() - t0)["elapsed"]))
    }

  }

  # ---- Collect results ----
  I_true_matrix <- matrix(I_true_vec, nrow = length(times), ncol = n_sites)
  results <- list(
    true_incidence         = I_true_matrix,
    observed_cases         = observed_cases,
    mechanistic_prediction = I_star,
    x_star                 = x_star,
    epsilon_true           = epsilon_true,
    coords                 = coords,
    prev_data              = prev_data,
    draws_gp_with          = draws_gp_with,
    draws_gp_without       = draws_gp_without,
    spatial_coords         = spatial_coords,
    times                  = times,
    true_params            = true_params
  )

  # ---- Diagnostic Plots ----
  site_idx <- 1

  # Plot 1: Stage 1 mechanistic prediction vs true (with GP residuals) vs observed
  plot_data <- data.frame(
    month            = seq_along(times),
    true_incidence   = I_true_matrix[, site_idx],
    observed_cases   = observed_cases[, site_idx],
    mechanistic_pred = I_star[, site_idx]
  )
  p1 <- ggplot(plot_data, aes(x = month)) +
    geom_line(aes(y = true_incidence,   colour = "True Incidence (with GP)"), linewidth = 1) +
    geom_point(aes(y = observed_cases,  colour = "Observed Cases"), alpha = 0.6, size = 2) +
    geom_line(aes(y = mechanistic_pred, colour = "Mechanistic I* (no GP)"),
              linetype = "dashed", linewidth = 1) +
    scale_colour_manual(values = c("True Incidence (with GP)" = "#2E75B6",
                                   "Observed Cases" = "black",
                                   "Mechanistic I* (no GP)" = "#C00000")) +
    labs(title = sprintf("Stage 1 vs True Incidence (Site %d)", site_idx),
         subtitle = "Gap between dashed (I*) and solid (true) = GP residual structure",
         x = "Month", y = "Monthly Incidence / Cases", colour = "") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"))
  print(p1)

  # Plot 2: True GP residual surface
  eps_matrix <- matrix(epsilon_true, nrow = length(times), ncol = n_sites)
  eps_df <- expand.grid(month = seq_along(times), site = seq_len(n_sites))
  eps_df$epsilon <- as.vector(eps_matrix)
  p2 <- ggplot(eps_df, aes(x = month, y = factor(site), fill = epsilon)) +
    geom_tile() +
    scale_fill_gradient2(low = "#C00000", mid = "white", high = "#2E75B6", midpoint = 0) +
    labs(title = "True GP Residuals (epsilon)", x = "Month", y = "Site", fill = "epsilon") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  print(p2)

  # Plot 3: Alpha-gamma joint posterior (if GP draws available)
  if (!is.null(draws_gp_with)) {
    gp_df <- as.data.frame(as.matrix(draws_gp_with))
    p3 <- ggplot(gp_df, aes(x = alpha, y = gamma_rr)) +
      geom_point(alpha = 0.15, size = 0.8, colour = "#2E75B6") +
      geom_vline(xintercept = true_params$alpha, linetype = "dashed", colour = "grey30") +
      geom_hline(yintercept = true_params$reporting_rate, linetype = "dashed", colour = "grey30") +
      labs(title = "Alpha-Gamma Joint Posterior (GP+offset)",
           subtitle = "Cluster = identifiable; ridge = non-identifiable",
           x = "alpha (intercept)", y = "gamma (reporting rate)") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(face = "bold"))
    print(p3)

    # Plot 4: GP hyperparameter posteriors
    gp_hyper_df <- data.frame(
      log_sigma = gp_df$log_sigma,
      log_phi   = gp_df$log_phi,
      log_theta = gp_df$log_theta
    )
    true_vals <- data.frame(
      param = c("log_sigma", "log_phi", "log_theta"),
      value = c(log(true_params$gp_sigma2), log(true_params$gp_phi), log(true_params$gp_theta))
    )

    gp_long <- tidyr::pivot_longer(gp_hyper_df, everything(), names_to = "param", values_to = "value")
    p4 <- ggplot(gp_long, aes(x = value)) +
      geom_density(fill = "#2E75B6", alpha = 0.4) +
      geom_vline(data = true_vals, aes(xintercept = value),
                 linetype = "dashed", colour = "#C00000") +
      facet_wrap(~ param, scales = "free") +
      labs(title = "GP Hyperparameter Posteriors", x = "Value", y = "Density") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(face = "bold"), strip.text = element_text(face = "bold"))
    print(p4)

    # Plot 5: Posterior predictive check
    gamma_est <- median(gp_df$gamma_rr)
    alpha_est <- median(gp_df$alpha)
    pred_gp_with <- gamma_est * exp(alpha_est) * I_star[, site_idx]

    pp_df <- data.frame(
      month    = seq_along(times),
      observed = observed_cases[, site_idx],
      gp_with  = pred_gp_with
    )

    p5 <- ggplot(pp_df, aes(x = month)) +
      geom_point(aes(y = observed, colour = "Observed"), size = 2, alpha = 0.7) +
      geom_line(aes(y = gp_with, colour = "GP+offset"), linewidth = 1) +
      scale_colour_manual(values = c("Observed" = "black", "GP+offset" = "#2E75B6")) +
      labs(title = sprintf("Posterior Predictive Check (Site %d)", site_idx),
           x = "Month", y = "Predicted Cases", colour = "") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
    print(p5)
  }

  return(results)
}


# ==============================================================================
# POSTERIOR EXTRACTION AND METRICS
# ==============================================================================

#' Extract posterior summary from MCMC draws.
#' @param draws MCMC draw object from greta::mcmc()
#' @param prob_lower,prob_upper Credible interval bounds (default 0.025, 0.975)
#' @return Data frame with mean, median, sd, lower, upper per parameter
#' @export
extract_posterior_summary <- function(draws, prob_lower = 0.025, prob_upper = 0.975) {
  draws_mat <- if (inherits(draws, "mcmc.list")) as.matrix(draws) else as.matrix(draws)
  draws_df <- as.data.frame(draws_mat)
  data.frame(
    parameter = colnames(draws_df),
    mean   = colMeans(draws_df, na.rm = TRUE),
    median = apply(draws_df, 2, median, na.rm = TRUE),
    sd     = apply(draws_df, 2, sd, na.rm = TRUE),
    lower  = apply(draws_df, 2, quantile, probs = prob_lower, na.rm = TRUE),
    upper  = apply(draws_df, 2, quantile, probs = prob_upper, na.rm = TRUE),
    row.names = NULL
  )
}


#' Compute performance metrics (RMSE, MAE, coverage) for model comparison.
#' @param predicted,truth Predicted and true values (vectors or matrices)
#' @param lower_ci,upper_ci Optional CI bounds for coverage calculation
#' @return List with rmse, mae, relative_error, coverage, n_obs
#' @export
compute_performance_metrics <- function(predicted, truth,
                                        lower_ci = NULL, upper_ci = NULL) {
  pred_vec  <- as.vector(predicted)
  truth_vec <- as.vector(truth)
  valid     <- !is.na(pred_vec) & !is.na(truth_vec)
  pred_vec  <- pred_vec[valid]
  truth_vec <- truth_vec[valid]

  rmse      <- sqrt(mean((pred_vec - truth_vec)^2))
  mae       <- mean(abs(pred_vec - truth_vec))
  rel_error <- mean(abs(pred_vec - truth_vec) / pmax(abs(truth_vec), 1e-6))

  coverage <- NA
  if (!is.null(lower_ci) && !is.null(upper_ci)) {
    lower_vec <- as.vector(lower_ci)[valid]
    upper_vec <- as.vector(upper_ci)[valid]
    coverage  <- mean(truth_vec >= lower_vec & truth_vec <= upper_vec)
  }

  list(rmse = rmse, mae = mae, relative_error = rel_error,
       coverage = coverage, n_obs = length(truth_vec))
}



# ==============================================================================
# MAIN EXAMPLE
# ==============================================================================

#' Complete workflow demonstration — GP + dual likelihood pipeline
#' @param n_sites Number of sites (default 10)
#' @param n_times Number of monthly time points (default 48)
#' @param n_samples,warmup,chains MCMC settings
#' @return Simulation results (invisibly)
#' @export
main_example <- function(n_sites = 10, n_times = 48,
                         n_samples = 1000, warmup = 500, chains = 2) {
  message("EpiWave FOI Model: GP + Dual Likelihood Demonstration")

  results <- simulate_and_estimate(
    n_sites               = n_sites,
    n_times               = n_times,
    include_interventions = TRUE,
    run_mcmc              = TRUE,
    n_samples             = n_samples,
    warmup                = warmup,
    chains                = chains,
    use_sparse_gp         = TRUE,
    n_inducing            = 25
  )

  return(invisible(results))
}


# ==============================================================================
# INTERACTIVE EXECUTION
# ==============================================================================

if (interactive()) {
  results <- main_example()
}
