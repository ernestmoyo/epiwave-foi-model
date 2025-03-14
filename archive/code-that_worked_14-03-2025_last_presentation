# Libraries
library(greta)
library(greta.dynamics)

# Ross-McDonald equations integrated clearly with Bayesian framework

# Differential equations:
dx_dt <- function(x, z, m, a, b, r) {
  m * a * b * z * (1 - x) - r * x
}

dz_dt <- function(x, z, a, c, g) {
  a * c * x * (1 - z) - g * z
}

# State iteration function
iterate_state <- function(state, iter, m, a, b, c, r, g) {
  x <- state[1]  # human infection proportion
  z <- state[2]  # mosquito infection proportion

  state_new <- c(
    x + dx_dt(x, z, m, a, b, r),
    z + dz_dt(x, z, a, c, g)
  )
  state_new
}

# function to simulate the model
simulate_integrated_malaria_model <- function(n_days = 365, n_sims = 100, X_covariate) {

  # Priors
  alpha <- normal(0, 1)                    # Intercept for mosquito abundance
  beta <- normal(1, 1, truncation = c(0, Inf))  # Covariate effect (positive constraint)

  # Mosquito abundance spatial/temporal model
  log_m <- alpha + beta * X_covariate
  m <- exp(log_m)

  # Additional parameter priors
  a <- normal(0.2, 0.05, truncation = c(0,1))  # human biting rate
  b <- normal(0.6, 0.05, truncation = c(0,1))  # infection prob. mosquito to human
  c <- normal(0.9, 0.05, truncation = c(0,1))  # infection prob. human to mosquito
  inv_r <- normal(14, 1, truncation = c(0, Inf))  # duration infection (humans)
  inv_g <- normal(7, 1, truncation = c(0, Inf))   # lifespan mosquitoes

  # Convert to rates
  r <- 1/inv_r
  g <- 1/inv_g

  # Initial conditions (small infection)
  initial_state <- c(0.01, 0)

  # Solve differential equations dynamically
  solution <- iterate_dynamic_function(
    transition_function = iterate_state,
    initial_state = initial_state,
    niter = n_days,
    tol = 0,
    m = m,
    a = a,
    b = b,
    c = c,
    r = r,
    g = g,
    parameter_is_time_varying = "m"
  )

  # Extract states clearly
  x <- solution$all_states[1, ]  # Humans infected
  z <- solution$all_states[2, ]  # Mosquitoes infected

  # Force of infection (FOI)
  log_S <- log(m * a * b * t(z))

  # Gaussian Process to capture unexplained spatial-temporal variability
  epsilon <- greta.gp::gp(matrix(1:n_days), kernel = greta.gp::rbf(lengthscale = 30,
                                                                   variance = 1))

  log_xi <- log_S + epsilon
  xi <- exp(log_xi)

  # Human infection incidence
  M_human <- normal(1000, 100)  # Example human population per location
  I_lt <- M_human * xi
  gamma_lt <- normal(0.1, 0.01) # example symptomatic rate

  # Poisson distributed cases
  C_lt <- poisson(I_lt * gamma_lt)

  # Model run: capture uncertainty clearly explained
  sims <- calculate(x, z, xi, C_lt, nsim = n_sims)

  # Output neatly explained plots
  par(mfrow=c(4,1))

  # Proportion of humans infected
  plot(sims$x[1,1,], type='n', ylim=c(0,1), main="Humans Infected (x)", ylab="Proportion infected", xlab="Day")
  for(i in 1:n_sims) lines(sims$x[i,1,])

  # Proportion of mosquitoes infected
  plot(sims$z[1,1,], type='n', ylim=c(0,1), main="Mosquitoes Infected (z)", ylab="Proportion infected", xlab="Day")
  for(i in 1:n_sims) lines(sims$z[i,1,])

  # Force of infection (xi)
  plot(sims$xi[1,1,], type='n', main="Force of Infection (xi)", ylab="FOI", xlab="Day")
  for(i in 1:n_sims) lines(sims$xi[i,1,])

  # Simulated Cases
  plot(sims$C_lt[1,1,], type='n', main="Simulated Malaria Cases (C_lt)", ylab="Cases", xlab="Day")
  for(i in 1:n_sims) lines(sims$C_lt[i,1,])

  # Return values clearly explained
  return(list(
    humans_infected = sims$x,
    mosquitoes_infected = sims$z,
    force_of_infection = sims$xi,
    simulated_cases = sims$C_lt
  ))
}

# Example covariate data (random here, replace with real covariates)
X_covariate <- sin(2*pi*(1:365)/365)

# usage
results <- simulate_integrated_malaria_model(n_days=365, n_sims=10, X_covariate=X_covariate)

