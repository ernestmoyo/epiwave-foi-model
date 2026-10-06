# ==============================================================================
# EpiWave FOI model: vector-informed malaria incidence mapping
# Ernest Moyo (NM-AIST / Vector Atlas), PhD Objective 2
# ==============================================================================
#
# The model has two stages:
#
#   Stage 1 (fixed, solved once per site)
#     m, a, g (with ITN/IRS effects) -> Ross-Macdonald ODE -> z
#     I* = m * a * b * z                          mechanistic incidence RATE
#
#   Stage 2 (Bayesian, greta)
#     log I = alpha + log I* + epsilon            I* = 0 drops the offset
#     epsilon_t = theta * epsilon_{t-1} + f_t,    f_t ~ GP(0, sigma2 * Matern52(phi))
#     cases      ~ Poisson(gamma * I * N)          N = population
#     positives  ~ Binomial(T, p),  p = 1 - exp(-(I convolved with q))
#                                                  T = number tested
#
# Every matrix in this file is [n_sites x n_times]: rows are sites, columns are
# months. This is the epiwave.mapping convention.
#
# Sections
#   1. Priors and constants
#   2. Stage 1: entomology, interventions, ODE, I*
#   3. Stage 2: AR(1), detectability kernel, the greta model
#   4. Simulating data from known truth
#   5. Scoring a fit against the truth
#
# The worked example (simulate, fit with and without I*, plot) is R/run_demo.R.
# The multi-replicate study is R/sim_estimation_harness.R.
# ==============================================================================

library(deSolve)
library(dplyr)
library(tidyr)


# ==============================================================================
# 1. PRIORS AND CONSTANTS
# ==============================================================================

# All Stage 2 priors, in one place. fit_epiwave_gp() and
# R/diagnostics/prior_predictive.R both read this list.
#
# The marginal variance tau2 is sampled and the innovation variance is derived,
# sigma2 = tau2 * (1 - theta^2). The data only see the marginal variance, so
# sampling sigma2 and theta directly gave a ridge the chains could not cross.
# See docs/2026-10-05_code_simplification/prior_notes.md.
#
# What each prior implies (95% prior range, on a scale that can be checked):
#   alpha     incidence is 0.14 to 7.1 times I* before residuals. On real data,
#             where I* from these toy inputs is far too high, this needs revisiting.
#   gamma_rr  a reporting rate of 0.017 to 0.20.
#   tau2      0.23 to 1.62, so a site one sd above the mechanism has 1.6 to 3.6
#             times I*.
#   theta     0.09 to 0.91: from little to strong month-to-month persistence.
#   phi       0.62 to 4.39 on the unit square: correlation between sites one unit
#             apart of 0.24 to 0.96.
PRIORS <- list(
  alpha    = list(mean = 0, sd = 1),                     # intercept
  gamma_rr = list(mean = 0.1, sd = 0.05, lower = 0.001), # reporting rate
  tau2     = list(meanlog = -0.5, sdlog = 0.5),          # stationary marginal variance of epsilon
  theta    = list(shape1 = 2, shape2 = 2),               # AR(1) correlation
  phi      = list(meanlog = 0.5, sdlog = 0.5)            # spatial lengthscale
)

I_STAR_FLOOR <- 1e-6   # keeps log(I*) finite where the ODE gives z ~ 0
DAYS_PER_STEP <- 30    # one model time step ("month") = 30 days, for the ODE and the kernel
MAX_DETECT_DAYS <- 30  # how long an infection stays detectable by the test


# ==============================================================================
# 2. STAGE 1: FIXED ENTOMOLOGY -> ROSS-MACDONALD -> I*
# ==============================================================================
# m, a and g are synthetic stand-ins for the simulation study. When the Vector
# Atlas surfaces arrive (posterior medians per pixel and month) they replace
# these three functions; nothing downstream changes. b, c and r are held
# fixed everywhere, since they are not expected to vary much in space or time.

# Mosquito-to-human ratio m, with two rainy seasons a year (East Africa).
# The same seasonal curve is used at every site.
get_fixed_m <- function(time, location, baseline_m = 2.0,
                        seasonal_amplitude = 0.5, phase_shift = 0) {
  year <- time / 365
  seasonal <- 1 + seasonal_amplitude *
    (sin(4 * pi * year + phase_shift) - 0.3 * cos(2 * pi * year))
  m_t <- baseline_m * pmax(seasonal, 0.1)
  matrix(m_t, nrow = length(location), ncol = length(time), byrow = TRUE)
}

# Human biting rate a (bites per mosquito per day). Vector Atlas gives m * a
# together; a is separated by dividing by a mapped a.
get_fixed_a <- function(time, location, baseline_a = 0.3) {
  matrix(baseline_a, nrow = length(location), ncol = length(time))
}

# m from the interim Vector Atlas abundance maps. Those maps are calibrated
# against human biting data, so they give m * a (bites per person per day);
# dividing by a mapped biting rate a gives m. The maps are observed under the
# nets already in use: applying ITN effects to this m again counts them twice.
m_from_biting_rate <- function(biting_rate, a) {
  stopifnot(all(a > 0))
  biting_rate / a
}

# Mosquito death rate g (per day), i.e. a 10-day mean lifespan.
get_fixed_g <- function(time, location, baseline_g = 1/10) {
  matrix(baseline_g, nrow = length(location), ncol = length(time))
}

# Share of the ITN effect kept when vectors are resistant, given the fraction
# susceptible in bioassays (the Vector Atlas IR-cube quantity). Posterior mean of
# the Symons et al. relationship, as used in goldingn/ir_cube
# (R/fig_ento_epi_impact.R): 100% coverage with a fully resistant population
# works like 54% coverage with a susceptible one, because nets still block bites.
# It was estimated on the prevalence scale, so applying it to coverage is an
# interim step until a hut-trial mapping (bioassay -> killing, deterrence) is used.
itn_effect_retained <- function(susceptibility) {
  1 - 0.46 * (1 - susceptibility)
}

# ITN and IRS act on m, a and g before the ODE, so interventions enter
# through the entomology rather than the statistical model. Each effect is a
# multiplier that equals 1 at zero coverage. Effect sizes follow Griffin et al.
# (2010) and Bhatt et al. (2015). Resistance lowers the effective coverage.
# IRS has its own susceptibility, because IRS insecticides are mostly not the
# pyrethroids used on nets; it scales the IRS effect linearly.
apply_interventions <- function(m, a, g,
                                itn_coverage = NULL, irs_coverage = NULL,
                                itn_susceptibility = 1, irs_susceptibility = 1,
                                itn_kill_rate = 0.5, itn_feeding_inhibit = 0.3,
                                itn_mortality_boost = 0.3,
                                irs_efficacy = 0.5, irs_feeding_inhibit = 0.1) {
  if (!is.null(itn_coverage)) {
    stopifnot(all(itn_coverage >= 0 & itn_coverage <= 1))
    n <- itn_coverage * itn_effect_retained(itn_susceptibility)
    m <- m * (1 - n * itn_kill_rate)
    a <- a * (1 - n * itn_feeding_inhibit)
    g <- g * (1 + n * itn_mortality_boost)
  }
  if (!is.null(irs_coverage)) {
    stopifnot(all(irs_coverage >= 0 & irs_coverage <= 1))
    s <- irs_coverage * irs_susceptibility
    g <- g * (1 + s * irs_efficacy)
    a <- a * (1 - s * irs_feeding_inhibit)
  }
  list(m = m, a = a, g = g)
}

# Extrinsic incubation period (days) from temperature: the degree-day form of
# Gething et al. (2011), as in modd-africa/hackthon2026 (EIP.R). Below 16 C the
# parasite cannot complete sporogony.
eip_gething <- function(temperature) {
  ifelse(temperature <= 16, Inf, 111 / (temperature - 16))
}

# Fraction of infected mosquitoes that survive the EIP when it is split into
# `stages` exposed stages: (1 + g n / stages)^-stages, which tends to exp(-g n).
# One stage gives 0.50 at g n = 1 against the exact 0.37; four stages give 0.41.
eip_survival <- function(g, eip_days, stages = 4) {
  (1 + g * eip_days / stages)^(-stages)
}

# Ross-Macdonald: x = human prevalence, z = mosquito infection prevalence.
# m, a and g arrive as functions of time so they can vary through the year.
ross_macdonald_ode <- function(t, state, parms) {
  x <- state[1]
  z <- state[2]
  m <- parms$m(t); a <- parms$a(t); g <- parms$g(t)

  dx_dt <- m * a * parms$b * z * (1 - x) - parms$r * x
  dz_dt <- a * parms$c * x * (1 - z) - g * z
  list(c(dx_dt, dz_dt))
}

# The same model with an extrinsic incubation period: newly infected mosquitoes
# pass through `stages` exposed classes y (each left at rate stages / n) before
# becoming infectious. Only (1 + g n / stages)^-stages of them live that long.
ross_macdonald_eip_ode <- function(t, state, parms) {
  k <- parms$stages
  x <- state[1]
  y <- state[1 + seq_len(k)]
  z <- state[k + 2]
  m <- parms$m(t); a <- parms$a(t); g <- parms$g(t)
  leave <- k * parms$inv_eip(t)   # 1 / EIP is 0 when it is too cold for sporogony

  dx_dt <- m * a * parms$b * z * (1 - x) - parms$r * x
  dy_dt <- c(a * parms$c * x * (1 - sum(y) - z), leave * y[-k]) - (g + leave) * y
  dz_dt <- leave * y[k] - g * z
  list(c(dx_dt, dy_dt, dz_dt))
}

# Solve the ODE once per site (no inference here: this is the point of the
# two-stage design). Returns x and z as [n_sites x n_times].
#   eip_days = NULL    no incubation period (the 2-state model above)
#   eip_days = number or [n_sites x n_times] matrix, e.g. eip_gething(temperature)
solve_ross_macdonald_multi_site <- function(m_matrix, a_matrix, g_matrix, times,
                                            b = 0.8, c = 0.8, r = 1/7,
                                            x0 = 0.01, z0 = 0.001,
                                            eip_days = NULL, eip_stages = 4) {
  stopifnot(identical(dim(m_matrix), dim(a_matrix)),
            identical(dim(m_matrix), dim(g_matrix)),
            ncol(m_matrix) == length(times))
  x <- z <- matrix(NA_real_, nrow = nrow(m_matrix), ncol = length(times))
  if (!is.null(eip_days)) {
    stopifnot(length(eip_days) %in% c(1, length(m_matrix)))
    eip_days <- matrix(eip_days, nrow = nrow(m_matrix), ncol = length(times))
  }

  for (site in seq_len(nrow(m_matrix))) {
    parms <- list(m = approxfun(times, m_matrix[site, ], rule = 2),
                  a = approxfun(times, a_matrix[site, ], rule = 2),
                  g = approxfun(times, g_matrix[site, ], rule = 2),
                  b = b, c = c, r = r)
    if (is.null(eip_days)) {
      solution <- ode(y = c(x = x0, z = z0), times = times,
                      func = ross_macdonald_ode, parms = parms, method = "lsoda")
    } else {
      eip <- eip_days[site, ]
      parms$inv_eip <- approxfun(times, 1 / eip, rule = 2)
      parms$stages <- eip_stages
      # start the exposed classes near balance with z0: about g n z0 in total
      y0 <- if (is.finite(eip[1])) g_matrix[site, 1] * eip[1] * z0 / eip_stages else 0
      state0 <- c(x = x0, setNames(rep(y0, eip_stages), paste0("y", seq_len(eip_stages))), z = z0)
      solution <- ode(y = state0, times = times,
                      func = ross_macdonald_eip_ode, parms = parms, method = "lsoda")
    }
    x[site, ] <- solution[, "x"]
    z[site, ] <- solution[, "z"]
  }
  list(x = x, z = z)
}

# I* = m * a * b * z, new human infections per person per day. It is a rate:
# population enters later, in the case likelihood.
compute_mechanistic_prediction <- function(m_matrix, a_matrix, b, z_matrix) {
  m_matrix * a_matrix * b * z_matrix
}


# ==============================================================================
# 3. STAGE 2: GP + AR(1) RESIDUALS, DUAL LIKELIHOOD
# ==============================================================================

# AR(1) through time, written as one matrix multiply so TensorFlow runs it fast:
# epsilon_t = sum_i rho^i * innovation_{t-i}.
# Taken unchanged from epiwave.mapping (R/ar1.R).
ar1 <- function(rho, innovations) {
  n_times <- ncol(innovations)
  t_seq <- seq_len(n_times)
  t_mat <- outer(t_seq, t_seq, FUN = "-")
  t_mat <- pmax(t_mat, 0)
  mask <- lower.tri(t_mat, diag = TRUE)
  rho_mat <- (rho ^ t_mat) * mask
  t(rho_mat %*% t(innovations))
}

# q(d): probability that a person infected d days ago still tests positive.
# A rise-then-decay placeholder, the same one epiwave.mapping's sim_data.R uses.
# Replace it with a published diagnostic-sensitivity curve for real data.
default_q_daily <- function(days, max_detect_days = MAX_DETECT_DAYS) {
  up   <- plogis(days * 2 - 5)
  down <- 1 - plogis(days / 2 - 10)
  q <- up * down
  q[days < 0 | days > max_detect_days] <- 0
  q
}

# Turn a daily kernel into weights per model time step (here, per month).
# From epiwave.mapping (R/transform_convolution_kernel.R), with two
# changes: the step length is rounded to whole days, and lags outside the
# kernel return 0 instead of NA.
transform_convolution_kernel <- function(kernel_daily, max_diff_days, timeperiod_days) {
  timeperiod_days <- round(timeperiod_days)
  max_timeperiods <- 1 + ceiling(max_diff_days / timeperiod_days)
  total_integral <- sum(kernel_daily(0:max_diff_days))

  table <- tidyr::expand_grid(
    infection_day = seq_len(max_timeperiods * timeperiod_days),
    test_day = seq_len(timeperiod_days)
  ) %>%
    mutate(day_diff = infection_day - test_day) %>%
    filter(day_diff >= 0) %>%
    mutate(kernel_val_day = kernel_daily(day_diff),
           timeperiod_diff = (infection_day - 1) %/% timeperiod_days) %>%
    group_by(test_day, timeperiod_diff) %>%
    summarise(partial_integral = sum(kernel_val_day), .groups = "drop") %>%
    mutate(fraction = partial_integral / total_integral) %>%
    group_by(timeperiod_diff) %>%
    summarise(fraction = mean(fraction), .groups = "drop") %>%
    mutate(kernel_val = fraction * total_integral)

  function(time_period_difference) {
    result <- table$kernel_val[match(time_period_difference, table$timeperiod_diff)]
    result[is.na(result)] <- 0
    result
  }
}

# The convolution as a fixed matrix C, so greta computes it as one multiply:
# (I %*% C)[, t] = sum_k I[, t - k] * weights[k + 1].
build_convolution_matrix <- function(n_times, weights) {
  lag <- outer(seq_len(n_times), seq_len(n_times), function(i, j) j - i)
  C <- matrix(0, nrow = n_times, ncol = n_times)
  in_kernel <- lag >= 0 & lag < length(weights)
  C[in_kernel] <- weights[lag[in_kernel] + 1]
  C
}

# The monthly detectability matrix used by both the simulator and the model.
build_detectability_matrix <- function(n_times) {
  q_monthly <- transform_convolution_kernel(kernel_daily = default_q_daily,
                                            max_diff_days = MAX_DETECT_DAYS,
                                            timeperiod_days = DAYS_PER_STEP)
  max_lag <- 1 + ceiling(MAX_DETECT_DAYS / DAYS_PER_STEP)
  build_convolution_matrix(n_times, q_monthly(0:max_lag))
}

# The Stage 2 model. All data matrices are [n_sites x n_times].
#
#   use_mechanistic = FALSE  drops log(I*): the standard geostatistical model
#                            (the I* = 0 comparison).
#   center_alpha = TRUE      subtracts the realised mean of epsilon, so alpha is
#                            the only intercept. This changes what alpha means:
#                            it then estimates alpha + mean(epsilon).
#   inducing                 optional [m x 2] inducing points for a sparse GP.
#   observed_sites           rows whose data enter the likelihood (NULL = all).
#                            The GP still covers every site, so the fit predicts
#                            the others: this is how held-out validation works.
#
# Prevalence data are required. Without them alpha and gamma are not
# separately identifiable, so there is deliberately no case-only path.
#
# Returns a greta model. attr(model, "I_latent") holds the latent incidence so
# a fit can be turned into a map with greta::calculate().
fit_epiwave_gp <- function(observed_cases, I_star, N_pop, spatial_coords,
                           prev_data, prev_conv_matrix,
                           use_mechanistic = TRUE, center_alpha = FALSE,
                           inducing = NULL, gp_tol = 1e-3,
                           observed_sites = NULL) {
  if (!"package:greta.gp" %in% search())
    stop("greta is not loaded. Run source('R/greta_setup.R') first.")
  stopifnot(identical(dim(observed_cases), dim(I_star)),
            identical(dim(observed_cases), dim(N_pop)))
  n_times <- ncol(observed_cases)

  # Priors (section 1). Declaration order fixes greta's parameter layout.
  alpha    <- normal(PRIORS$alpha$mean, PRIORS$alpha$sd)
  gamma_rr <- normal(PRIORS$gamma_rr$mean, PRIORS$gamma_rr$sd,
                     truncation = c(PRIORS$gamma_rr$lower, Inf))
  tau2     <- lognormal(PRIORS$tau2$meanlog, PRIORS$tau2$sdlog)
  theta    <- beta(PRIORS$theta$shape1, PRIORS$theta$shape2)
  sigma2   <- tau2 * (1 - theta ^ 2)
  phi      <- lognormal(PRIORS$phi$meanlog, PRIORS$phi$sdlog)

  # Residuals: a spatial GP draw for each month, chained through time by AR(1).
  kernel  <- mat52(lengthscales = phi, variance = sigma2)
  f       <- gp(spatial_coords, kernel, inducing = inducing, n = n_times, tol = gp_tol)
  epsilon <- ar1(rho = theta, innovations = f)
  if (center_alpha) epsilon <- epsilon - mean(epsilon)

  # Latent infection incidence (per person per day)
  log_I <- if (use_mechanistic) {
    alpha + log(pmax(I_star, I_STAR_FLOOR)) + epsilon
  } else {
    alpha + epsilon
  }
  I_latent <- exp(log_I)

  # Cases: population enters here, not in I*
  expected_cases <- gamma_rr * I_latent * N_pop
  if (is.null(observed_sites)) {
    cases <- as_data(observed_cases)
    distribution(cases) <- poisson(expected_cases)
  } else {
    cases <- as_data(observed_cases[observed_sites, , drop = FALSE])
    distribution(cases) <- poisson(expected_cases[observed_sites, ])
  }

  # Prevalence depends on I, not on the ODE's x, so the GP residuals and alpha
  # inform both likelihoods. Recent infections are summed through the
  # detectability kernel q. epiwave.mapping uses the linear sum p = sum(I * q);
  # here p = 1 - exp(-sum(I * q)), the chance of at least one detectable
  # infection. The two agree when the sum is small, and this form cannot
  # exceed 1.
  survey_used <- if (is.null(observed_sites)) seq_along(prev_data$survey_indices) else
    which(((prev_data$survey_indices - 1) %% nrow(observed_cases) + 1) %in% observed_sites)
  stopifnot(length(survey_used) > 0)   # the surveys identify alpha and gamma
  detectable <- I_latent %*% as_data(prev_conv_matrix)
  p_positive <- 1 - exp(-detectable[prev_data$survey_indices[survey_used]])
  positives  <- as_data(prev_data$n_positive[survey_used])
  distribution(positives) <- binomial(prev_data$n_tested[survey_used], p_positive)

  greta_model <- model(alpha, gamma_rr, sigma2, phi, theta)
  attr(greta_model, "I_latent") <- I_latent
  greta_model
}


# ==============================================================================
# 4. SIMULATING DATA FROM KNOWN TRUTH
# ==============================================================================

# The truth the simulation study tries to recover.
# gp_phi = 3 keeps the model sampleable. The prior predictive check shows phi
# itself is not identifiable at this value (see the paper).
TRUE_PARAMS <- list(
  baseline_m = 2.0, baseline_a = 0.3, baseline_g = 1/10,
  b = 0.8, c = 0.8, r = 1/7,
  population = 10000, reporting_rate = 0.1,
  itn_max = 0.7,              # ITN coverage reached by the final month
  itn_susceptibility = 0.8,   # fraction susceptible in bioassays (IR-cube scale)
  alpha = 0, gp_sigma = 0.6, gp_phi = 3.0, gp_rho = 0.75
)

# True residuals: Matern 5/2 in space, AR(1) in time. Built with plain R
# (MASS::mvrnorm) rather than greta, so the truth does not share code with the
# model being tested.
simulate_gp_residuals <- function(spatial_coords, n_times, sigma, phi, rho) {
  n_sites <- nrow(spatial_coords)
  r <- sqrt(5) * as.matrix(dist(spatial_coords)) / phi
  K_space <- (1 + r + r^2 / 3) * exp(-r) + diag(1e-6, n_sites)

  innovations <- matrix(NA_real_, nrow = n_sites, ncol = n_times)
  for (t in seq_len(n_times)) {
    innovations[, t] <- MASS::mvrnorm(1, mu = rep(0, n_sites), Sigma = K_space)
  }

  epsilon <- matrix(0, nrow = n_sites, ncol = n_times)
  epsilon[, 1] <- sigma * innovations[, 1]
  for (t in 2:n_times) {
    epsilon[, t] <- rho * epsilon[, t - 1] + sigma * innovations[, t]
  }
  epsilon
}

# Prevalence surveys at a random share of site-months. survey_indices index
# as.vector() of the [n_sites x n_times] prevalence matrix.
simulate_prevalence_surveys <- function(prevalence, survey_fraction = 0.3,
                                        sample_size_range = c(50, 200),
                                        seed = 789) {
  set.seed(seed)
  n_surveys <- round(length(prevalence) * survey_fraction)
  survey_indices <- sort(sample(seq_along(prevalence), n_surveys))
  true_prev <- pmax(as.vector(prevalence)[survey_indices], 1e-6)
  n_tested <- sample(sample_size_range[1]:sample_size_range[2], n_surveys, replace = TRUE)

  list(survey_indices = survey_indices,
       n_tested = n_tested,
       n_positive = rbinom(n_surveys, size = n_tested, prob = true_prev),
       true_prevalence = true_prev)
}

# One full dataset: Stage 1, true residuals, cases and surveys.
# seed = NULL reproduces the demo; the harness passes one seed per replicate.
# n_times is the number of monthly steps; the matrices have n_times + 1
# columns because month 0 is included.
# site_variation = TRUE gives each site its own ITN coverage (east-west) and
# susceptibility (south-north), so I* varies in space as it will with real
# intervention and IR-cube maps. FALSE gives every site the same entomology.
simulate_epiwave_data <- function(n_sites = 10, n_times = 48,
                                  true_params = TRUE_PARAMS,
                                  include_interventions = TRUE, seed = NULL,
                                  site_variation = FALSE) {
  tp <- true_params
  seeds <- if (is.null(seed)) c(coords = 123, eps = 321, cases = 456, surveys = 789)
           else seed + c(coords = 0L, eps = 1L, cases = 2L, surveys = 3L)

  times     <- seq(0, n_times * DAYS_PER_STEP, by = DAYS_PER_STEP)
  locations <- sprintf("Site_%02d", seq_len(n_sites))

  set.seed(seeds[["coords"]])
  coords <- matrix(runif(n_sites * 2, min = -5, max = 5), ncol = 2,
                   dimnames = list(locations, c("lon", "lat")))
  coords_01 <- apply(coords, 2, function(v) (v - min(v)) / (diff(range(v)) + 1e-10))

  # Stage 1
  m <- get_fixed_m(times, locations, baseline_m = tp$baseline_m, seasonal_amplitude = 0.6)
  a <- get_fixed_a(times, locations, baseline_a = tp$baseline_a)
  g <- get_fixed_g(times, locations, baseline_g = tp$baseline_g)
  if (include_interventions) {
    itn_max <- if (site_variation) 0.2 + 0.7 * coords_01[, "lon"] else rep(tp$itn_max, n_sites)
    susceptibility <- if (site_variation) 0.4 + 0.6 * coords_01[, "lat"] else tp$itn_susceptibility
    itn_coverage <- outer(itn_max, seq(0, 1, length.out = length(times)))
    adjusted <- apply_interventions(
      m, a, g,
      itn_coverage = itn_coverage,
      itn_susceptibility = matrix(susceptibility, nrow = n_sites, ncol = length(times))
    )
    m <- adjusted$m; a <- adjusted$a; g <- adjusted$g
  }
  ode <- solve_ross_macdonald_multi_site(m, a, g, times = times,
                                         b = tp$b, c = tp$c, r = tp$r)
  I_star <- compute_mechanistic_prediction(m, a, tp$b, ode$z)

  # True incidence = I* adjusted by the true residuals
  set.seed(seeds[["eps"]])
  epsilon <- simulate_gp_residuals(coords_01, n_times = length(times),
                                   sigma = tp$gp_sigma, phi = tp$gp_phi, rho = tp$gp_rho)
  I_true <- exp(tp$alpha + log(pmax(I_star, I_STAR_FLOOR)) + epsilon)

  # Observations
  population <- matrix(tp$population, nrow = n_sites, ncol = length(times))
  set.seed(seeds[["cases"]])
  expected_cases <- tp$reporting_rate * I_true * population
  observed_cases <- matrix(rpois(length(expected_cases), expected_cases),
                           nrow = n_sites, ncol = length(times))

  conv_matrix <- build_detectability_matrix(length(times))
  prevalence_true <- 1 - exp(-(I_true %*% conv_matrix))
  prev_data <- simulate_prevalence_surveys(prevalence_true, seed = seeds[["surveys"]])

  list(observed_cases = observed_cases, I_star = I_star, x_star = ode$x,
       pop_matrix = population, spatial_coords_norm = coords_01,
       conv_matrix = conv_matrix, prev_data = prev_data,
       I_true_mat = I_true, prevalence_true_mat = prevalence_true,
       epsilon_true_mat = epsilon, times = times, true_params = tp)
}


# ==============================================================================
# 5. SCORING A FIT AGAINST THE TRUTH
# ==============================================================================

# Posterior mean, median, sd and 95% interval for each parameter.
extract_posterior_summary <- function(draws, prob_lower = 0.025, prob_upper = 0.975) {
  d <- as.data.frame(as.matrix(draws))
  data.frame(parameter = colnames(d),
             mean   = colMeans(d),
             median = apply(d, 2, median),
             sd     = apply(d, 2, sd),
             lower  = apply(d, 2, quantile, probs = prob_lower),
             upper  = apply(d, 2, quantile, probs = prob_upper),
             row.names = NULL)
}

# How close a predicted surface is to the truth. If interval bounds are given,
# coverage is the share of true values inside them.
compute_performance_metrics <- function(predicted, truth, lower_ci = NULL, upper_ci = NULL) {
  error <- as.vector(predicted) - as.vector(truth)
  coverage <- if (is.null(lower_ci)) NA else
    mean(truth >= lower_ci & truth <= upper_ci)
  list(rmse = sqrt(mean(error^2)),
       mae = mean(abs(error)),
       relative_error = mean(abs(error) / pmax(abs(as.vector(truth)), 1e-6)),
       coverage = coverage,
       n_obs = length(error))
}
