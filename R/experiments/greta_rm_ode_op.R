
# =============================================================================
# GRETA OPERATION: Ross-MacDonald ODE Solver
# =============================================================================
# Wraps the TensorFlow DormandPrince ODE solver as a greta operation
# for use in Bayesian inference models.
#
# Author: Ernest Moyo (building on Nick Golding's integrate_RMd)
# Date: 2025-01-24
# =============================================================================

library(greta)
library(tensorflow)

# Source the base functions if not already loaded
if (!exists("integrate_RMd")) {
  source("R/experiments/tf_ode_test.R")
}

#' Create a greta operation for Ross-MacDonald ODE integration
#'
#' This function creates greta arrays for the ODE solution, allowing
#' Bayesian inference on model parameters (m_int, m_slope, etc.)
#'
#' @param m_int Intercept for log(m) adjustment (greta array)
#' @param m_slope Slope for log(m) adjustment (greta array)
#' @param a_int Intercept for logit(a) adjustment (greta array)
#' @param a_slope Slope for logit(a) adjustment (greta array)
#' @param g_int Intercept for log(g) adjustment (greta array)
#' @param g_slope Slope for log(g) adjustment (greta array)
#' @param x_0 Initial human infection proportion (scalar)
#' @param z_0 Initial mosquito infection proportion (scalar)
#' @param b Transmission probability mosquito->human (fixed scalar)
#' @param c Transmission probability human->mosquito (fixed scalar)
#' @param r Human recovery rate (fixed scalar)
#' @param knots Spline knot locations (vector)
#' @param m_star_coefs Pre-computed spline coefficients for m (matrix)
#' @param a_star_coefs Pre-computed spline coefficients for a (matrix)
#' @param g_star_coefs Pre-computed spline coefficients for g (matrix)
#' @param burnin Burn-in days before evaluation period
#' @param times Times at which to evaluate the solution
#' @param rtol Relative tolerance for ODE solver
#' @param atol Absolute tolerance for ODE solver
#'
#' @return A greta array with dimensions [n_times, n_sites, 2] containing
#'         x (human prevalence) and z (mosquito prevalence) over time
#'
greta_ross_macdonald <- function(
    # Greta arrays for inference parameters
    m_int, m_slope,
    a_int, a_slope, 
    g_int, g_slope,
    # Fixed parameters
    x_0 = 0.1, z_0 = 0.1,
    b = 0.8, c = 0.8, r = 1/7,
    # Pre-computed spline data
    knots, m_star_coefs, a_star_coefs, g_star_coefs,
    # Time parameters
    burnin = 365, times = seq(0, 4*365, by = 30),
    # Solver parameters
    rtol = 0.001, atol = 1e-06
) {
  
  # Get dimensions
  n_sites <- ncol(m_star_coefs)
  n_times <- length(times)
  
  # Define the tensorflow operation
  tf_op <- function(m_int, m_slope, a_int, a_slope, g_int, g_slope) {
    
    # Call the integrate_RMd function
    integral <- integrate_RMd(
      burnin = burnin,
      times = times,
      x_0 = x_0, z_0 = z_0,
      m_int = m_int, m_slope = m_slope,
      a_int = a_int, a_slope = a_slope,
      g_int = g_int, g_slope = g_slope,
      b = b, c = c, r = r,
      knots = knots,
      m_star_coefs = m_star_coefs,
      a_star_coefs = a_star_coefs,
      g_star_coefs = g_star_coefs,
      rtol = rtol, atol = atol
    )
    
    integral
  }
  
  # Create the greta operation
  op <- greta::.internals$nodes$constructors$op
  
  result <- op(
    tf_op,
    m_int, m_slope, a_int, a_slope, g_int, g_slope,
    operation_args = list(),
    tf_operation = "custom_rm_ode",
    dim = c(n_times, n_sites, 2L)
  )
  
  result
}


#' Extract human prevalence (x) from ODE solution
#' @param solution Output from greta_ross_macdonald
#' @return Greta array of human prevalence [n_times, n_sites]
get_human_prevalence <- function(solution) {
  solution[, , 1]
}

#' Extract mosquito prevalence (z) from ODE solution
#' @param solution Output from greta_ross_macdonald
#' @return Greta array of mosquito prevalence [n_times, n_sites]
get_mosquito_prevalence <- function(solution) {
  solution[, , 2]
}

#' Calculate Force of Infection from ODE solution
#' @param solution Output from greta_ross_macdonald
#' @param m_t Mosquito abundance at evaluation times (matrix: n_times x n_sites)
#' @param a_t Biting rate at evaluation times (matrix: n_times x n_sites)
#' @param b Transmission probability mosquito->human
#' @return Greta array of FOI [n_times, n_sites]
get_foi <- function(solution, m_t, a_t, b = 0.8) {
  z <- get_mosquito_prevalence(solution)
  m_t * a_t * b * z
}


# =============================================================================
# EXAMPLE USAGE
# =============================================================================

example_greta_model <- function() {
  
  cat("Setting up example greta model...\n")
  
  # --- Pre-compute spline data (done once, outside the model) ---
  
  # Simulation settings
  days_evaluate <- seq(0, 4*365, by = 30)
  burnin_days <- 365
  n_sites <- 5
  
  # Setup observation times
  previous_days_observed <- -rev(seq(0, burnin_days, by = 30)[-1])
  days_observed <- c(previous_days_observed, days_evaluate)
  
  # Mock observed data (in practice, from Vector Atlas / entomological surveys)
  m_star_observed <- t(replicate(n_sites, m(days_observed)))
  a_star_observed <- t(replicate(n_sites, a(days_observed)))
  g_star_observed <- t(replicate(n_sites, g(days_observed)))
  
  # Define knots
  n_knots <- length(days_observed)
  epsilon <- 1e-3
  limits <- c(-burnin_days - epsilon, max(days_evaluate) + epsilon)
  knots <- define_knots(n_knots, limits)
  
  # Pre-compute spline coefficients
  m_star_coefs <- get_spline_coefs(log(m_star_observed), days_observed, knots)
  a_star_coefs <- get_spline_coefs(qlogis(a_star_observed), days_observed, knots)
  g_star_coefs <- get_spline_coefs(log(g_star_observed), days_observed, knots)
  
  # --- Define greta model ---
  
  # Priors on adjustment parameters
  m_int <- normal(0, 0.5)
  m_slope <- normal(1, 0.2, truncation = c(0, Inf))
  a_int <- normal(0, 0.5)
  a_slope <- normal(1, 0.2, truncation = c(0, Inf))
  g_int <- normal(0, 0.5)
  g_slope <- normal(1, 0.2, truncation = c(0, Inf))
  
  # Solve ODEs
  solution <- greta_ross_macdonald(
    m_int = m_int, m_slope = m_slope,
    a_int = a_int, a_slope = a_slope,
    g_int = g_int, g_slope = g_slope,
    knots = knots,
    m_star_coefs = m_star_coefs,
    a_star_coefs = a_star_coefs,
    g_star_coefs = g_star_coefs,
    burnin = burnin_days,
    times = days_evaluate
  )
  
  # Extract prevalence
  x_pred <- get_human_prevalence(solution)
  
  # --- Connect to observed data ---
  # (In practice, link to malaria case data via observation model)
  # observed_cases ~ poisson(population * incidence_rate)
  
  cat("Model setup complete!\n")
  cat("Dimensions of solution:", dim(solution), "\n")
  
  # Return model components
  list(
    solution = solution,
    x_pred = x_pred,
    params = list(m_int = m_int, m_slope = m_slope,
                  a_int = a_int, a_slope = a_slope,
                  g_int = g_int, g_slope = g_slope)
  )
}

cat("Greta operation wrapper loaded!\n")
cat("Use greta_ross_macdonald() to create ODE solution in greta models.\n")
cat("Run example_greta_model() for a demonstration.\n")

