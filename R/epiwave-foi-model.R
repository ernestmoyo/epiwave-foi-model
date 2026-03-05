# ==============================================================================
# EpiWave FOI Model: Vector-Informed Malaria Transmission Mapping
# ==============================================================================
# Author: Ernest Moyo (NM-AIST / Vector Atlas)
# PhD Objective 2 — Framework Development
#
# Two-stage computational strategy:
#   Stage 1: Fixed entomological params (m,a,g) -> Ross-Macdonald ODE -> I*(s,t)
#   Stage 2: NegBin likelihood with log(I*) as offset, 2 MCMC params (log_rate, log_size)
#
# Key design choice (per supervisor Prof Nick Golding):
#   Do NOT infer dynamic ODE parameters. Use fixed Vector Atlas / temperature
#   model estimates, solve ODEs ONCE per pixel, then calibrate to case data.
# ==============================================================================

suppressPackageStartupMessages({
  library(deSolve)
  library(ggplot2)
  library(tidyverse)
})
# greta loaded lazily inside fit_epiwave_with_offset()

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
# STAGE 2: BAYESIAN CALIBRATION (Negative Binomial with mechanistic offset)
# ==============================================================================

#' Fit Stage 2 NegBin model with/without mechanistic offset via greta.
#'
#' Model: log(mu) = log_rate + log(I*), C ~ NegBin(size, prob)
#' Priors: log_rate ~ N(-2,1), log_size ~ N(3,1), size = exp(log_size)
#'
#' @param observed_cases Matrix [n_times x n_sites] of case counts
#' @param I_star Matrix [n_times x n_sites] of Stage 1 mechanistic predictions
#' @param use_mechanistic If TRUE, uses log(I*) as offset; if FALSE, intercept-only baseline
#' @return greta model object for mcmc()
#' @export
fit_epiwave_with_offset <- function(observed_cases, I_star,
                                    use_mechanistic = TRUE) {
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
  suppressPackageStartupMessages(library(greta))

  cases_vec  <- as.vector(observed_cases)
  I_star_vec <- as.vector(I_star)

  # Priors (only 2 free parameters)
  log_rate <- normal(-2, 1)
  log_size <- normal(3, 1)
  size     <- exp(log_size)

  # Linear predictor
  log_mu <- if (use_mechanistic) {
    log_rate + log(I_star_vec + 1e-10)
  } else {
    log_rate + log(mean(I_star_vec) + 1e-10)
  }

  # Likelihood
  mu   <- exp(log_mu)
  prob <- size / (size + mu)
  distribution(cases_vec) <- negative_binomial(size, prob)

  model(log_rate, log_size)
}


# ==============================================================================
# VALIDATION: SIMULATION-ESTIMATION STUDY
# ==============================================================================

#' Run simulation-estimation study comparing WITH vs WITHOUT mechanistic offset.
#' @param n_sites Number of spatial locations (default 10)
#' @param n_times Number of monthly time points (default 48)
#' @param true_params List of true parameter values (NULL for defaults)
#' @param include_interventions Simulate ITN scale-up (default TRUE)
#' @param run_mcmc Whether to run MCMC sampling (default TRUE)
#' @param n_samples,warmup,chains MCMC settings
#' @param draws_with,draws_without Pre-computed MCMC draws (optional)
#' @return List of simulation results
#' @export
simulate_and_estimate <- function(n_sites = 10, n_times = 48,
                                  true_params = NULL,
                                  include_interventions = TRUE,
                                  run_mcmc = TRUE,
                                  n_samples = 1000, warmup = 500, chains = 2,
                                  draws_with = NULL, draws_without = NULL) {

  # Step 1: Generate synthetic data
  if (is.null(true_params)) {
    true_params <- list(
      baseline_m = 2.0, baseline_a = 0.3, baseline_g = 1/10,
      b = 0.8, c = 0.8, r = 1/7,
      population = 10000, reporting_rate = 0.1
    )
  }

  times     <- seq(0, n_times * 30, by = 30)
  locations <- paste0("Site_", sprintf("%02d", 1:n_sites))

  set.seed(123)
  spatial_coords <- matrix(runif(n_sites * 2, min = -5, max = 5),
                           ncol = 2, dimnames = list(locations, c("lon", "lat")))

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

  I_true <- compute_mechanistic_prediction(
    m_matrix = m_true, a_matrix = a_true, b = true_params$b,
    z_matrix = ode_solution$z, population_matrix = pop_matrix
  )

  expected_cases <- true_params$reporting_rate * I_true
  set.seed(456)
  observed_cases <- matrix(rpois(length(expected_cases), lambda = expected_cases),
                           nrow = length(times), ncol = n_sites)

  I_star_mech <- compute_mechanistic_prediction(
    m_matrix = m_true, a_matrix = a_true, b = true_params$b,
    z_matrix = ode_solution$z, population_matrix = pop_matrix
  )

  # Step 2: Build and sample WITH offset model
  model_with_mech <- fit_epiwave_with_offset(
    observed_cases = observed_cases, I_star = I_star_mech, use_mechanistic = TRUE
  )
  if (run_mcmc && is.null(draws_with)) {
    t0 <- proc.time()
    draws_with <- tryCatch(
      greta::mcmc(model_with_mech, n_samples = n_samples,
                  warmup = warmup, chains = chains, verbose = FALSE),
      error = function(e) { message("MCMC error (WITH): ", e$message); NULL }
    )
    message(sprintf("WITH model: %.0fs", (proc.time() - t0)["elapsed"]))
  }

  # Step 3: Build and sample WITHOUT offset model
  model_without_mech <- fit_epiwave_with_offset(
    observed_cases = observed_cases, I_star = I_star_mech, use_mechanistic = FALSE
  )
  if (run_mcmc && is.null(draws_without)) {
    t0 <- proc.time()
    draws_without <- tryCatch(
      greta::mcmc(model_without_mech, n_samples = n_samples,
                  warmup = warmup, chains = chains, verbose = FALSE),
      error = function(e) { message("MCMC error (WITHOUT): ", e$message); NULL }
    )
    message(sprintf("WITHOUT model: %.0fs", (proc.time() - t0)["elapsed"]))
  }

  # Step 4: Collect results
  results <- list(
    true_incidence         = I_true,
    observed_cases         = observed_cases,
    mechanistic_prediction = I_star_mech,
    model_with_mech        = model_with_mech,
    model_without_mech     = model_without_mech,
    draws_with             = draws_with,
    draws_without          = draws_without,
    spatial_coords         = spatial_coords,
    times                  = times,
    true_params            = true_params
  )

  # Step 5: Plots
  site_idx <- 1

  # Plot 1: Stage 1 mechanistic prediction vs observed
  plot_data <- data.frame(
    month            = seq_along(times),
    true_incidence   = I_true[, site_idx],
    observed_cases   = observed_cases[, site_idx],
    mechanistic_pred = I_star_mech[, site_idx]
  )
  p1 <- ggplot(plot_data, aes(x = month)) +
    geom_line(aes(y = true_incidence,   colour = "True Incidence"), linewidth = 1) +
    geom_point(aes(y = observed_cases,  colour = "Observed Cases"), alpha = 0.6, size = 2) +
    geom_line(aes(y = mechanistic_pred, colour = "Mechanistic Prediction I*"),
              linetype = "dashed", linewidth = 1) +
    scale_colour_manual(values = c("True Incidence" = "#2E75B6",
                                   "Observed Cases" = "black",
                                   "Mechanistic Prediction I*" = "#C00000")) +
    labs(title = sprintf("Stage 1 -- Mechanistic Prediction vs Observed (Site %d)", site_idx),
         x = "Month", y = "Monthly Incidence / Cases", colour = "") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"))
  print(p1)

  # Plots 2-5: Stage 2 posterior (only if MCMC draws available)
  if (!is.null(draws_with) && !is.null(draws_without)) {
    draws_with_df    <- as.data.frame(as.matrix(draws_with))
    draws_without_df <- as.data.frame(as.matrix(draws_without))
    draws_with_df$model    <- "WITH mechanistic offset"
    draws_without_df$model <- "WITHOUT mechanistic offset"
    draws_combined <- rbind(draws_with_df, draws_without_df)

    colour_vals <- c("WITH mechanistic offset" = "#2E75B6",
                     "WITHOUT mechanistic offset" = "#C00000")

    # Plot 2: Posterior density -- log_rate
    p2 <- ggplot(draws_combined, aes(x = log_rate, fill = model, colour = model)) +
      geom_density(alpha = 0.35, linewidth = 0.8) +
      geom_vline(xintercept = log(0.1), linetype = "dashed", colour = "grey30") +
      scale_fill_manual(values = colour_vals) +
      scale_colour_manual(values = colour_vals) +
      labs(title = "Posterior: log_rate", x = "log_rate", y = "Density", fill = "", colour = "") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
    print(p2)

    # Plot 3: Posterior density -- log_size
    p3 <- ggplot(draws_combined, aes(x = log_size, fill = model, colour = model)) +
      geom_density(alpha = 0.35, linewidth = 0.8) +
      scale_fill_manual(values = colour_vals) +
      scale_colour_manual(values = colour_vals) +
      labs(title = "Posterior: log_size (overdispersion)",
           x = "log_size", y = "Density", fill = "", colour = "") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
    print(p3)

    # Plot 4: MCMC trace plots
    n_iter <- nrow(draws_with_df) / 2
    trace_df <- rbind(
      data.frame(iteration = rep(seq_len(n_iter), 2),
                 chain = rep(c("Chain 1", "Chain 2"), each = n_iter),
                 log_rate = as.matrix(draws_with)[, "log_rate"],
                 log_size = as.matrix(draws_with)[, "log_size"],
                 model = "WITH mechanistic offset"),
      data.frame(iteration = rep(seq_len(n_iter), 2),
                 chain = rep(c("Chain 1", "Chain 2"), each = n_iter),
                 log_rate = as.matrix(draws_without)[, "log_rate"],
                 log_size = as.matrix(draws_without)[, "log_size"],
                 model = "WITHOUT mechanistic offset")
    )

    for (param in c("log_rate", "log_size")) {
      p4 <- ggplot(trace_df, aes(x = iteration, y = .data[[param]],
                                  colour = chain, alpha = chain)) +
        geom_line(linewidth = 0.4) +
        scale_alpha_manual(values = c("Chain 1" = 0.9, "Chain 2" = 0.6)) +
        scale_colour_manual(values = c("Chain 1" = "#2E75B6", "Chain 2" = "#C00000")) +
        facet_wrap(~ model, ncol = 1) +
        labs(title = paste("MCMC Trace --", param), x = "Iteration", y = param) +
        theme_minimal(base_size = 11) +
        theme(legend.position = "bottom", strip.text = element_text(face = "bold"),
              plot.title = element_text(face = "bold"))
      print(p4)
    }

    # Plot 5: Posterior predictive check
    reporting_rate_with    <- exp(median(draws_with_df$log_rate))
    reporting_rate_without <- exp(median(draws_without_df$log_rate))
    pp_df <- data.frame(
      month        = seq_along(times),
      observed     = observed_cases[, site_idx],
      pred_with    = reporting_rate_with    * I_star_mech[, site_idx],
      pred_without = reporting_rate_without * mean(I_star_mech)
    )
    p5 <- ggplot(pp_df, aes(x = month)) +
      geom_point(aes(y = observed, colour = "Observed Cases"), size = 2, alpha = 0.7) +
      geom_line(aes(y = pred_with, colour = "WITH offset"), linewidth = 1) +
      geom_line(aes(y = pred_without, colour = "WITHOUT offset"),
                linewidth = 1, linetype = "dashed") +
      scale_colour_manual(values = c("Observed Cases" = "black",
                                     "WITH offset" = "#2E75B6",
                                     "WITHOUT offset" = "#C00000")) +
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
  draws_df <- if (is.data.frame(draws)) draws else as.data.frame(draws)
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

#' Complete workflow demonstration — Stage 1 + Stage 2 plots from cache
#' @return Simulation results (invisibly)
#' @export
main_example <- function(n_samples = 1000, warmup = 500, chains = 2) {


  cat("==============================================================================\n")
  cat("  EpiWave FOI Model: Complete Workflow Demonstration\n")
  cat("==============================================================================\n\n")


  results <- simulate_and_estimate(
    n_sites               = 10,
    n_times               = 48,
    include_interventions = TRUE,
    run_mcmc              = TRUE,
    n_samples             = n_samples,
    warmup                = warmup,
    chains                = chains
  )

  return(invisible(results))
}


# ==============================================================================
# INTERACTIVE EXECUTION
# ==============================================================================

# Runs automatically when sourced. MCMC takes ~ 10 seconds. All 5 plots on completion.
if (interactive()) {
  results <- main_example()
}
