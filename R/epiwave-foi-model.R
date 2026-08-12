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
#   - Use a GP for spatially-correlated residuals.
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


#' Compute mechanistic infection incidence rate: I*(t,s) = m*a*b*z
#'
#' Per Nick Golding's specification: I* is a RATE, not a count.
#' Population enters the Poisson likelihood separately.
#'
#' @param m_matrix,a_matrix,z_matrix Matrices [n_times x n_sites]
#' @param b Transmission probability mosquito->human
#' @return Matrix I_star [n_times x n_sites] — infection incidence rate
#' @export
compute_mechanistic_prediction <- function(m_matrix, a_matrix, b, z_matrix) {
  if (!all(dim(m_matrix) == dim(a_matrix),
           dim(m_matrix) == dim(z_matrix)))
    stop("All input matrices must have identical dimensions")

  m_matrix * a_matrix * b * z_matrix
}


# ==============================================================================
# STAGE 2: BAYESIAN CALIBRATION (GP residuals + dual likelihood)
# ==============================================================================
# Model specification (per Nick Golding / epiwave.mapping):
#   log(I_{l,t}) = alpha + log(I*_{l,t}) + epsilon_{l,t}
#   epsilon ~ GP(0, K), K = sigma^2 * K_space(phi) * K_time(theta)
#   C_{l,t} ~ Poisson(gamma * I_{l,t} * N_pop_{l})     -- case counts (N_pop = population)
#   Y_{l,t} ~ Binomial(n_test_{l,t}, prev_{l,t})       -- prevalence surveys (n_test = sample size)
# Prevalence depends on the GP-adjusted incidence I (NOT mechanistic X), via the
# epiwave.mapping detectability convolution:
#   intensity_{l,t} = sum_k I_{l,t-k} * w_k   (w from a test-detectability kernel)
#   prev_{l,t}      = 1 - exp(-intensity_{l,t})
# This makes the prevalence likelihood depend on alpha, so the dual likelihood
# can separate alpha (mean incidence) from gamma (reporting rate).
# ==============================================================================

#' Build spatial-only GP kernel (Matern 5/2).
#'
#' Temporal correlation is handled separately via AR(1), following epiwave.mapping.
#'
#' @param phi Spatial lengthscale (greta scalar)
#' @param sigma2 GP marginal variance (greta scalar)
#' @return greta.gp kernel object (spatial only)
#' @export
build_gp_kernel <- function(phi, sigma2) {
  if (!requireNamespace("greta.gp", quietly = TRUE))
    stop("Package 'greta.gp' required. Install with install.packages('greta.gp')")
  greta.gp::mat52(lengthscales = phi, variance = sigma2)
}


#' AR(1) temporal correlation for GP innovations.
#'
#' Ported from epiwave.mapping/R/ar1.R (Nick Golding).
#' Expands the iterative AR(1) equation X_t = rho * X_{t-1} + eps_t
#' into a matrix multiplication for efficient TensorFlow execution.
#'
#' @param rho AR(1) correlation coefficient (greta scalar, 0 to 1)
#' @param innovations Matrix [n_sites x n_times] of spatial GP innovations
#' @return Matrix [n_sites x n_times] of temporally correlated residuals
#' @keywords internal
ar1 <- function(rho, innovations) {
  n_times <- ncol(innovations)
  t_seq <- seq_len(n_times)
  t_mat <- outer(t_seq, t_seq, FUN = "-")
  t_mat <- pmax(t_mat, 0)
  mask <- lower.tri(t_mat, diag = TRUE)
  rho_mat <- (rho ^ t_mat) * mask
  t(rho_mat %*% t(innovations))
}


#' Daily test-detectability kernel for the infection -> prevalence convolution.
#'
#' q_daily(d) = probability an infected person still tests positive d days after
#' infection; zero outside [0, max_detect_days]. Default is a smooth
#' rise-then-decay curve (matching epiwave.mapping's sim_data.R placeholder) for
#' the simulation-estimation study. Swap for a literature diagnostic-sensitivity
#' curve when fitting real data.
#'
#' @param days Numeric vector of days since infection
#' @param max_detect_days Detectability horizon in days (default 30)
#' @return Numeric vector of detection probabilities in [0, 1]
#' @export
default_q_daily <- function(days, max_detect_days = 30) {
  up   <- plogis(days * 2 - 5)        # rises over the first few days
  down <- 1 - plogis(days / 2 - 10)   # decays toward the horizon
  q <- up * down
  q[days < 0] <- 0
  q[days > max_detect_days] <- 0
  q
}


#' Transform a daily convolution kernel to a discrete coarser-timestep kernel.
#'
#' Ported from epiwave.mapping/R/transform_convolution_kernel.R (Nick Golding).
#' Integrates a daily kernel into per-timeperiod weights, returning a function
#' that maps an integer timeperiod lag to its kernel weight.
#'
#' @param kernel_daily Function: vector of day-lags -> densities
#' @param max_diff_days Maximum day-lag with non-zero density
#' @param timeperiod_days Length of the modelling timestep in days
#' @return Function: integer timeperiod-lag -> kernel weight (0 outside support)
#' @export
transform_convolution_kernel <- function(kernel_daily, max_diff_days, timeperiod_days) {
  tp_days <- round(timeperiod_days)
  max_timeperiods <- 1 + ceiling(max_diff_days / tp_days)
  total_integral <- sum(kernel_daily(0:max_diff_days))

  table <- tidyr::expand_grid(
    infection_day = seq_len(max_timeperiods * tp_days),
    test_day      = seq_len(tp_days)
  ) %>%
    dplyr::mutate(day_diff = .data$infection_day - .data$test_day) %>%
    dplyr::filter(.data$day_diff >= 0) %>%
    dplyr::mutate(
      kernel_val_day  = kernel_daily(.data$day_diff),
      timeperiod_diff = (.data$infection_day - 1) %/% tp_days
    ) %>%
    dplyr::group_by(.data$test_day, .data$timeperiod_diff) %>%
    dplyr::summarise(partial_integral = sum(.data$kernel_val_day), .groups = "drop") %>%
    dplyr::mutate(fraction = .data$partial_integral / total_integral) %>%
    dplyr::group_by(.data$timeperiod_diff) %>%
    dplyr::summarise(fraction = mean(.data$fraction), .groups = "drop") %>%
    dplyr::mutate(kernel_val = .data$fraction * total_integral)

  function(time_period_difference) {
    idx <- match(time_period_difference, table$timeperiod_diff)
    result <- table$kernel_val[idx]
    result[is.na(result)] <- 0
    result
  }
}


#' Build a fixed [n_times x n_times] temporal convolution matrix.
#'
#' Lets the detectability convolution run as a single matrix multiply in greta
#' (the same trick ar1() uses): prevalence_intensity = I_latent %*% C reproduces
#'   intensity[, t] = sum_k I[, t - k] * weights[k + 1].
#'
#' @param n_times Number of time steps
#' @param weights Numeric vector; weights[k + 1] is the lag-k kernel weight
#' @return Numeric matrix [n_times x n_times], lower-banded (C[i, j] = w_{j - i})
#' @export
build_convolution_matrix <- function(n_times, weights) {
  K <- length(weights) - 1
  C <- matrix(0, nrow = n_times, ncol = n_times)
  for (lag in 0:K) {
    w <- weights[lag + 1]
    if (w == 0) next
    for (j in seq_len(n_times)) {
      i <- j - lag
      if (i >= 1) C[i, j] <- w
    }
  }
  C
}


#' Simulate GP residuals using spatial Matern 5/2 + AR(1) temporal.
#'
#' Matches the epiwave.mapping structure: spatial innovations from mvrnorm,
#' then AR(1) forward in time.
#'
#' @param spatial_coords Matrix [n_sites x 2] of (lon, lat)
#' @param n_times Number of time steps
#' @param sigma GP marginal SD (scalar)
#' @param phi Spatial lengthscale (scalar)
#' @param rho AR(1) temporal correlation (scalar, 0 to 1)
#' @return Matrix [n_sites x n_times] of GP residuals
#' @export
simulate_gp_residuals <- function(spatial_coords, n_times, sigma, phi, rho) {
  if (!requireNamespace("MASS", quietly = TRUE))
    stop("Package 'MASS' required for simulate_gp_residuals(). Install with install.packages('MASS')")
  n_sites <- nrow(spatial_coords)
  spatial_dists <- as.matrix(dist(spatial_coords))

  # Matern 5/2 spatial kernel (unit variance — sigma applied after AR(1))
  r_s <- sqrt(5) * spatial_dists / phi
  K_space <- (1 + r_s + r_s^2 / 3) * exp(-r_s) + diag(1e-6, n_sites)

  # Spatial innovations: one draw per time step
  innovations <- matrix(NA, nrow = n_sites, ncol = n_times)
  for (t in seq_len(n_times)) {
    innovations[, t] <- MASS::mvrnorm(1, mu = rep(0, n_sites), Sigma = K_space)
  }

  # Apply AR(1) and marginal variance
  epsilon <- matrix(0, nrow = n_sites, ncol = n_times)
  epsilon[, 1] <- sigma * innovations[, 1]
  for (t in 2:n_times) {
    epsilon[, t] <- rho * epsilon[, t - 1] + sigma * innovations[, t]
  }

  epsilon
}


#' Simulate Binomial prevalence survey data from a prevalence surface.
#'
#' Generates Binomial survey data at a random subset of (site, time) pairs. The
#' caller passes the I-derived, detectability-convolved prevalence surface
#' (1 - exp(-conv)), NOT mechanistic ODE prevalence X — matching the fitted
#' likelihood, which also depends on the GP-adjusted incidence I.
#'
#' @param x_matrix Matrix [n_sites x n_times] prevalence surface (I-derived)
#' @param gp_adjustment Vector of exp(epsilon) adjustments (length n_sites*n_times)
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
#' Follows epiwave.mapping: spatial GP innovations + AR(1) temporal correlation.
#'   log(I) = alpha + log(I*) + epsilon
#'   f ~ GP(0, sigma^2 * Matern52(phi)),  n = n_times spatial innovations
#'   epsilon = AR1(rho = theta, innovations = f)
#'   cases ~ Poisson(gamma * I * N)
#'   prevalence ~ Binomial(N_tested, 1 - exp(-(I convolved with detectability)))
#'
#' @param observed_cases Matrix [n_times x n_sites] of case counts
#' @param I_star Matrix [n_times x n_sites] of mechanistic incidence RATE
#' @param N_pop Matrix [n_times x n_sites] of population at risk (enters the
#'   Poisson case likelihood; distinct from n_test, the Binomial survey sample size)
#' @param spatial_coords Matrix [n_sites x 2] of (lon, lat) coordinates
#' @param prev_data List from simulate_prevalence_surveys() (NULL for case-only)
#' @param prev_conv_matrix Numeric [n_times x n_times] detectability convolution
#'   matrix from build_convolution_matrix(); required when prev_data is supplied
#' @param use_mechanistic If TRUE, use log(I*) as offset; if FALSE, drop offset (I*=0)
#' @param inducing Optional matrix [m x 2] of spatial inducing point coordinates
#' @param gp_tol Jitter for GP numerical stability (default 1e-3)
#' @param center_alpha If TRUE, subtract the realised mean of the GP+AR(1) field
#'   so alpha is the sole intercept (improves alpha mixing). Default FALSE.
#' @return greta model object for mcmc(); also carries attr "I_latent" (the
#'   latent incidence greta array, for greta::calculate posterior prediction)
#'   and attr "nodes" (the parameter greta arrays).
#' @export
fit_epiwave_gp <- function(observed_cases, I_star, N_pop,
                           spatial_coords, prev_data = NULL,
                           prev_conv_matrix = NULL,
                           use_mechanistic = TRUE,
                           inducing = NULL,
                           gp_tol = 1e-3,
                           center_alpha = FALSE) {
  if (!is.matrix(observed_cases))
    stop("observed_cases must be a matrix [n_times x n_sites]")
  if (!identical(dim(observed_cases), dim(I_star)))
    stop("observed_cases and I_star must have identical dimensions")

  if (!requireNamespace("greta", quietly = TRUE))
    stop("Package 'greta' required. Run R/greta_setup.R first.")
  if (!requireNamespace("greta.gp", quietly = TRUE))
    stop("Package 'greta.gp' required.")
  suppressPackageStartupMessages({
    library(greta)
    library(greta.gp)
  })

  n_times <- nrow(observed_cases)
  n_sites <- ncol(observed_cases)

  # --- Priors ---
  # NOTE (Aug 2026): making phi identifiable (prior median ~0.5, true 0.5) broke
  # MCMC mixing badly. The Workbench check traced this to a variance ridge between
  # sigma2 and theta (the AR(1) marginal variance is roughly sigma2 / (1 - theta^2),
  # so the two trade off), with theta having a flat prior that pins nothing. That is
  # a separate problem from phi. Reverted to the working config (phi median 1.65,
  # true 3.0). The next fix is to reparameterise that ridge and put a prior on theta.
  alpha     <- normal(0, 1)                                    # intercept
  gamma_rr  <- normal(0.1, 0.05, truncation = c(0.001, Inf))  # reporting rate
  sigma2    <- lognormal(-0.5, 0.5)                            # GP variance
  phi       <- lognormal(0.5, 0.5)                             # spatial lengthscale
  theta     <- variable(lower = 0, upper = 1)                  # AR(1) temporal correlation

  # --- Spatial GP + AR(1) temporal (epiwave.mapping pattern) ---
  K_space <- build_gp_kernel(phi, sigma2)
  coord_array <- as_data(spatial_coords)

  # Spatial innovations: n_times columns, each spatially correlated
  f <- if (!is.null(inducing)) {
    inducing_array <- as_data(inducing)
    gp(coord_array, K_space, inducing = inducing_array, n = n_times, tol = gp_tol)
  } else {
    gp(coord_array, K_space, n = n_times, tol = gp_tol)
  }

  # AR(1) temporal correlation: epsilon_mat is [n_sites x n_times]
  epsilon_mat <- ar1(rho = theta, innovations = f)

  # Optional centring: remove the realised field mean so alpha is the sole
  # intercept. Decouples alpha from the GP level and greatly improves alpha's
  # mixing/ESS. (alpha then estimates true_alpha + mean(epsilon); the harness
  # compares against that effective intercept.)
  if (center_alpha) epsilon_mat <- epsilon_mat - mean(epsilon_mat)

  # --- Latent infection incidence rate ---
  # I_star is a RATE [n_times x n_sites], transpose to [n_sites x n_times]
  log_I_mat <- if (use_mechanistic) {
    I_star_t <- t(pmax(I_star, 1e-6))
    alpha + log(I_star_t) + epsilon_mat
  } else {
    # I* = 0: standard geostatistical model, no mechanistic information
    alpha + epsilon_mat
  }
  I_latent_mat <- exp(log_I_mat)

  # --- Case likelihood (Poisson) ---
  # N_pop (population at risk) enters here, not in I* (Nick's specification)
  N_pop_t <- t(N_pop)
  cases_t <- t(observed_cases)
  cases_greta <- as_data(cases_t)
  expected_cases <- gamma_rr * I_latent_mat * N_pop_t
  distribution(cases_greta) <- poisson(expected_cases)

  # --- Prevalence likelihood (Binomial) if data provided ---
  # Prevalence depends on the GP-adjusted incidence I (NOT mechanistic X), via the
  # epiwave.mapping detectability convolution: intensity = I convolved with the
  # test-detectability kernel; prevalence = 1 - exp(-intensity) (bounded in (0,1)
  # and HMC-safe). Because I carries alpha + epsilon, this resolves the
  # alpha-gamma non-identifiability. Survey indices are into as.vector([n_sites x
  # n_times]) — the same column-major order as epsilon_mat / I_latent_mat.
  if (!is.null(prev_data)) {
    if (is.null(prev_conv_matrix))
      stop("prev_conv_matrix is required when prev_data is supplied")
    prev_intensity_mat <- I_latent_mat %*% as_data(prev_conv_matrix)  # [n_sites x n_times]
    intensity_vec <- prev_intensity_mat[seq_len(n_sites * n_times)]
    prev_prob     <- 1 - exp(-intensity_vec[prev_data$survey_indices])

    n_pos_greta <- as_data(prev_data$n_positive)
    distribution(n_pos_greta) <- binomial(prev_data$n_tested, prev_prob)
  }

  gm <- model(alpha, gamma_rr, sigma2, phi, theta)
  # Expose latent nodes (for posterior prediction via greta::calculate) without
  # changing the return type — greta::mcmc() ignores these attributes.
  attr(gm, "I_latent") <- I_latent_mat        # [n_sites x n_times], greta array
  attr(gm, "nodes") <- list(alpha = alpha, gamma_rr = gamma_rr,
                            sigma2 = sigma2, phi = phi, theta = theta)
  gm
}


# ==============================================================================
# VALIDATION: SIMULATION-ESTIMATION STUDY
# ==============================================================================

#' Simulate one EpiWave dataset (Stage 1 ODE + true GP residuals + observations).
#'
#' Extracted so the single-run demo and the multi-replicate harness share one
#' data-generating process. With seed = NULL the original fixed seeds are used
#' (reproduces the demo); pass a seed to get an independent replicate.
#'
#' @param n_sites,n_times Grid dimensions
#' @param true_params List of true values (NULL for defaults)
#' @param include_interventions Simulate ITN scale-up (default TRUE)
#' @param seed Integer base seed for this replicate (NULL = demo defaults)
#' @return List of everything needed to fit and to score: observed_cases, I_star,
#'   x_star, population, spatial_coords_norm, conv_matrix, prev_data, I_true_mat,
#'   prevalence_true_mat, epsilon_true_mat, times, true_params
#' @export
simulate_epiwave_data <- function(n_sites = 10, n_times = 48, true_params = NULL,
                                  include_interventions = TRUE, seed = NULL) {
  if (is.null(true_params)) {
    true_params <- list(
      baseline_m = 2.0, baseline_a = 0.3, baseline_g = 1/10,
      b = 0.8, c = 0.8, r = 1/7,
      population = 10000, reporting_rate = 0.1,
      # gp_phi = 3.0 is the working config: it keeps the model sampleable, so the
      # study recovers alpha and gamma (the main result). Making phi identifiable
      # (0.5) exposed a sigma2/theta variance ridge that broke mixing; reverted
      # pending a reparameterisation of that ridge.
      alpha = 0, gp_sigma = 0.6, gp_phi = 3.0, gp_rho = 0.75
    )
  }
  # Per-replicate seeds (NULL => original demo seeds)
  s_coord <- if (is.null(seed)) 123 else seed
  s_eps   <- if (is.null(seed)) 321 else seed + 1L
  s_case  <- if (is.null(seed)) 456 else seed + 2L
  s_surv  <- if (is.null(seed)) 789 else seed + 3L

  times     <- seq(0, n_times * 30, by = 30)
  locations <- paste0("Site_", sprintf("%02d", 1:n_sites))

  # Detectability convolution matrix (shared by truth and likelihood)
  month_days <- 365.25 / 12; max_detect_days <- 30
  q_monthly_fun <- transform_convolution_kernel(
    kernel_daily    = function(d) default_q_daily(d, max_detect_days = max_detect_days),
    max_diff_days   = max_detect_days, timeperiod_days = month_days)
  K_lag <- 1 + ceiling(max_detect_days / month_days)
  conv_matrix <- build_convolution_matrix(length(times), q_monthly_fun(0:K_lag))

  set.seed(s_coord)
  spatial_coords <- matrix(runif(n_sites * 2, min = -5, max = 5),
                           ncol = 2, dimnames = list(locations, c("lon", "lat")))

  m_true <- get_fixed_m(times, locations, baseline_m = true_params$baseline_m,
                        seasonal_amplitude = 0.6)
  a_true <- get_fixed_a(times, locations, baseline_a = true_params$baseline_a)
  g_true <- get_fixed_g(times, locations, baseline_g = true_params$baseline_g)

  if (include_interventions) {
    itn_coverage <- matrix(rep(seq(0, 0.7, length.out = length(times)), n_sites),
                           nrow = length(times), ncol = n_sites)
    adj <- apply_interventions(m = m_true, a = a_true, g = g_true,
                               itn_coverage = itn_coverage, resistance_index = 0.2)
    m_true <- adj$m; a_true <- adj$a; g_true <- adj$g
  }

  ode_solution <- solve_ross_macdonald_multi_site(
    m_matrix = m_true, a_matrix = a_true, g_matrix = g_true,
    times = times, b = true_params$b, c = true_params$c, r = true_params$r)

  pop_matrix <- matrix(true_params$population, nrow = length(times), ncol = n_sites)
  I_star <- compute_mechanistic_prediction(m_true, a_true, true_params$b, ode_solution$z)
  x_star <- ode_solution$x

  spatial_coords_norm <- cbind(
    lon = (spatial_coords[, 1] - min(spatial_coords[, 1])) /
      max(diff(range(spatial_coords[, 1])) + 1e-10),
    lat = (spatial_coords[, 2] - min(spatial_coords[, 2])) /
      max(diff(range(spatial_coords[, 2])) + 1e-10))

  set.seed(s_eps)
  epsilon_true_mat <- simulate_gp_residuals(
    spatial_coords = spatial_coords_norm, n_times = length(times),
    sigma = true_params$gp_sigma, phi = true_params$gp_phi, rho = true_params$gp_rho)

  I_star_t <- t(pmax(I_star, 1e-6)); pop_t <- t(pop_matrix)
  I_true_mat <- exp(true_params$alpha + log(I_star_t) + epsilon_true_mat)

  set.seed(s_case)
  expected_cases_mat <- true_params$reporting_rate * I_true_mat * pop_t
  cases_mat <- matrix(rpois(length(expected_cases_mat), expected_cases_mat),
                      nrow = n_sites, ncol = length(times))
  observed_cases <- t(cases_mat)

  prevalence_true_mat <- 1 - exp(-(I_true_mat %*% conv_matrix))
  prev_data <- simulate_prevalence_surveys(
    x_matrix = prevalence_true_mat, gp_adjustment = NULL,
    survey_fraction = 0.3, seed = s_surv)

  list(observed_cases = observed_cases, I_star = I_star, x_star = x_star,
       pop_matrix = pop_matrix, spatial_coords_norm = spatial_coords_norm,
       conv_matrix = conv_matrix, prev_data = prev_data,
       I_true_mat = I_true_mat, prevalence_true_mat = prevalence_true_mat,
       epsilon_true_mat = epsilon_true_mat, times = times, true_params = true_params)
}


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
#' @param seed Base seed for the simulated dataset (NULL = demo defaults)
#' @param center_alpha Centre the GP field so alpha is the sole intercept (default FALSE)
#' @return List of simulation results, including (when MCMC runs) posterior_summary
#'   and comparison_metrics (WITH offset vs I*=0 latent-incidence RMSE/MAE)
#' @export
simulate_and_estimate <- function(n_sites = 10, n_times = 48,
                                  true_params = NULL,
                                  include_interventions = TRUE,
                                  run_mcmc = TRUE,
                                  n_samples = 1000, warmup = 1000, chains = 2,
                                  use_sparse_gp = TRUE,
                                  n_inducing = 40,
                                  draws_gp_with = NULL,
                                  draws_gp_without = NULL,
                                  seed = NULL,
                                  center_alpha = FALSE) {

  # ---- Steps 1-6: simulate one dataset (shared with the harness) ----
  data <- simulate_epiwave_data(n_sites = n_sites, n_times = n_times,
                                true_params = true_params,
                                include_interventions = include_interventions,
                                seed = seed)
  true_params         <- data$true_params
  times               <- data$times
  spatial_coords_norm <- data$spatial_coords_norm
  conv_matrix         <- data$conv_matrix
  I_star              <- data$I_star
  x_star              <- data$x_star
  pop_matrix          <- data$pop_matrix
  epsilon_true_mat    <- data$epsilon_true_mat
  I_true_mat          <- data$I_true_mat
  prevalence_true_mat <- data$prevalence_true_mat
  observed_cases      <- data$observed_cases
  prev_data           <- data$prev_data

  message(sprintf("Data generated: %d sites x %d times, %d prevalence surveys",
                  n_sites, length(times), length(prev_data$n_positive)))

  # ---- Step 7: Spatial inducing points for sparse GP ----
  # Inducing points only add sparsity when there are fewer of them than sites.
  inducing <- NULL
  if (use_sparse_gp && n_inducing < n_sites) {
    ind_idx <- seq(1, n_sites, length.out = n_inducing)
    inducing <- spatial_coords_norm[round(ind_idx), , drop = FALSE]
    message(sprintf("Sparse GP: %d spatial inducing points", nrow(inducing)))
  } else if (use_sparse_gp) {
    message(sprintf(
      "Full GP: n_inducing (%d) >= n_sites (%d); inducing points add no sparsity, using full GP",
      n_inducing, n_sites))
  }

  # ---- Step 8: Fit GP+offset model ----
  if (run_mcmc) {
    message("\n--- Fitting GP + offset model ---")
    model_gp_with <- tryCatch(
      fit_epiwave_gp(
        observed_cases = observed_cases, I_star = I_star,
        N_pop = pop_matrix, spatial_coords = spatial_coords_norm,
        prev_data = prev_data, prev_conv_matrix = conv_matrix,
        use_mechanistic = TRUE, inducing = inducing, center_alpha = center_alpha
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
        observed_cases = observed_cases, I_star = I_star,
        N_pop = pop_matrix, spatial_coords = spatial_coords_norm,
        prev_data = prev_data, prev_conv_matrix = conv_matrix,
        use_mechanistic = FALSE, inducing = inducing, center_alpha = center_alpha
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
  # Transpose back to [n_times x n_sites] for storage/plotting
  I_true_matrix <- t(I_true_mat)
  results <- list(
    true_incidence         = I_true_matrix,
    observed_cases         = observed_cases,
    mechanistic_prediction = I_star,
    x_star                 = x_star,
    prevalence_true        = t(prevalence_true_mat),
    conv_matrix            = conv_matrix,
    epsilon_true_mat       = epsilon_true_mat,
    spatial_coords         = spatial_coords_norm,
    prev_data              = prev_data,
    draws_gp_with          = draws_gp_with,
    draws_gp_without       = draws_gp_without,
    times                  = times,
    true_params            = true_params
  )

  # ---- Diagnostic Plots ----
  if (!requireNamespace("gridExtra", quietly = TRUE))
    stop("Package 'gridExtra' required for the demo plots. Install with install.packages('gridExtra')")
  site_idx <- 1

  # Plot 1: Two-panel — rate space (I* vs I_true) and count space (expected vs observed)
  df_rate <- data.frame(
    month = seq_along(times),
    mechanistic_istar = I_star[, site_idx],
    true_incidence    = I_true_matrix[, site_idx]
  )
  df_count <- data.frame(
    month          = seq_along(times),
    expected_cases = true_params$reporting_rate * I_true_matrix[, site_idx] * pop_matrix[, site_idx],
    observed_cases = observed_cases[, site_idx]
  )
  p1a <- ggplot(df_rate, aes(x = month)) +
    geom_line(aes(y = mechanistic_istar, colour = "I* (mechanistic rate)"),
              linetype = "dashed", linewidth = 1) +
    geom_line(aes(y = true_incidence, colour = "I_true (with GP residuals)"), linewidth = 1) +
    scale_colour_manual(values = c("I* (mechanistic rate)" = "#C00000",
                                   "I_true (with GP residuals)" = "#2E75B6")) +
    labs(title = sprintf("Site %d — Rate Space", site_idx),
         subtitle = "Gap between lines = exp(epsilon) from spatial GP + AR(1)",
         x = NULL, y = "Infection incidence rate", colour = "") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"), axis.text.x = element_blank())
  p1b <- ggplot(df_count, aes(x = month)) +
    geom_line(aes(y = expected_cases, colour = "Expected (gamma * I_true * N)"), linewidth = 1) +
    geom_point(aes(y = observed_cases, colour = "Observed (Poisson draws)"),
               alpha = 0.6, size = 2) +
    scale_colour_manual(values = c("Expected (gamma * I_true * N)" = "#2E75B6",
                                   "Observed (Poisson draws)" = "black")) +
    labs(title = "Count Space",
         subtitle = "What Stage 2 actually fits to (Poisson likelihood)",
         x = "Month", y = "Case counts", colour = "") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"))
  gridExtra::grid.arrange(p1a, p1b, ncol = 1)

  # Plot 2: True GP residual surface
  # epsilon_true_mat is [n_sites x n_times], transpose for plotting
  eps_df <- expand.grid(month = seq_along(times), site = seq_len(n_sites))
  eps_df$epsilon <- as.vector(t(epsilon_true_mat))
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
      sigma2 = gp_df$sigma2,
      phi    = gp_df$phi,
      theta_rho = gp_df$theta
    )
    true_vals <- data.frame(
      param = c("sigma2", "phi", "theta_rho"),
      value = c(true_params$gp_sigma^2, true_params$gp_phi, true_params$gp_rho)
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

    # Plot 5: Posterior predictive check — predicted cases from the FITTED latent
    # incidence I (alpha + log(I*) + epsilon), via greta::calculate on attr
    # "I_latent", not just the mechanistic mean. Falls back to the mechanistic
    # mean only if the live model object is unavailable (e.g. draws passed in).
    gamma_est <- median(gp_df$gamma_rr)
    have_model_with <- exists("model_gp_with", inherits = FALSE) && !is.null(model_gp_with)
    if (have_model_with) {
      I_with_mat <- matrix(colMeans(as.matrix(greta::calculate(
        attr(model_gp_with, "I_latent"), values = draws_gp_with))),
        nrow = n_sites, ncol = length(times))            # [n_sites x n_times]
      pred_gp_with <- gamma_est * I_with_mat[site_idx, ] * pop_matrix[, site_idx]
    } else {
      alpha_est    <- median(gp_df$alpha)
      pred_gp_with <- gamma_est * exp(alpha_est) * I_star[, site_idx] * pop_matrix[, site_idx]
    }

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

    # Plot 6: Prevalence surveys — observed vs true at survey cells
    # Confirms the prevalence likelihood data is now derived from I (via the
    # detectability convolution), not from mechanistic X.
    prev_ppc_df <- data.frame(
      true_prevalence     = prev_data$true_prevalence,
      observed_prevalence = prev_data$n_positive / prev_data$n_tested
    )
    p6 <- ggplot(prev_ppc_df, aes(x = true_prevalence, y = observed_prevalence)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
      geom_point(alpha = 0.6, size = 2, colour = "#2E7D32") +
      labs(title = "Prevalence Surveys: Observed vs True",
           subtitle = "Prevalence derived from I via detectability convolution (1 - exp(-conv))",
           x = "True prevalence", y = "Observed (n_positive / n_tested)") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(face = "bold"))
    print(p6)

    # ---- Parameter recovery + WITH vs I*=0 comparison (Nick's headline test) ----
    param_summary <- extract_posterior_summary(draws_gp_with)
    message("\nGP + offset posterior summary:")
    print(param_summary, row.names = FALSE)
    results$posterior_summary <- param_summary

    if (have_model_with) {
      metrics_with <- compute_performance_metrics(t(I_with_mat), I_true_matrix)
      cmp <- data.frame(model = "GP + I* offset", rmse = metrics_with$rmse,
                        mae = metrics_with$mae, relative_error = metrics_with$relative_error)

      have_model_without <- exists("model_gp_without", inherits = FALSE) &&
        !is.null(model_gp_without) && !is.null(draws_gp_without)
      I_without_mat <- NULL
      if (have_model_without) {
        I_without_mat <- matrix(colMeans(as.matrix(greta::calculate(
          attr(model_gp_without, "I_latent"), values = draws_gp_without))),
          nrow = n_sites, ncol = length(times))
        metrics_without <- compute_performance_metrics(t(I_without_mat), I_true_matrix)
        cmp <- rbind(cmp, data.frame(model = "I* = 0 (geostatistical)",
                                     rmse = metrics_without$rmse, mae = metrics_without$mae,
                                     relative_error = metrics_without$relative_error))
      }
      message("\nLatent incidence recovery vs truth (lower RMSE = better):")
      print(cmp, row.names = FALSE)
      results$comparison_metrics <- cmp

      # Plot 7: predicted incidence at the example site — truth vs both models
      cmp_df <- data.frame(month = seq_along(times),
                           truth = I_true_matrix[, site_idx],
                           with_offset = I_with_mat[site_idx, ])
      cols <- c("Truth" = "black", "GP + I* offset" = "#2E75B6")
      p7 <- ggplot(cmp_df, aes(x = month)) +
        geom_line(aes(y = truth, colour = "Truth"), linewidth = 1) +
        geom_line(aes(y = with_offset, colour = "GP + I* offset"), linewidth = 1)
      if (have_model_without) {
        cmp_df$without_offset <- I_without_mat[site_idx, ]
        p7 <- p7 + geom_line(data = cmp_df,
                             aes(y = without_offset, colour = "I* = 0 (geostatistical)"),
                             linewidth = 1, linetype = "dashed")
        cols <- c(cols, "I* = 0 (geostatistical)" = "#C00000")
      }
      p7 <- p7 +
        scale_colour_manual(values = cols) +
        labs(title = sprintf("Latent Incidence Recovery (Site %d)", site_idx),
             subtitle = "WITH mechanistic offset vs I* = 0 (standard geostatistical)",
             x = "Month", y = "Infection incidence rate", colour = "") +
        theme_minimal(base_size = 12) +
        theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
      print(p7)
    }
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
  draws_df <- as.data.frame(as.matrix(draws))
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
                         n_samples = 2000, warmup = 1000, chains = 2) {
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
