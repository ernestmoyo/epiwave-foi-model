# Exports the real numbers the study site shows, plus reference values that the
# site's JavaScript re-implementation is tested against.
#   Rscript site/export_data.R      (from the project root)

source("R/epiwave-foi-model.R")
load("outputs/sim_estimation_results.RData")   # study

r4 <- function(x) signif(x, 6)

# Replicate 1 of the 50-replicate study (seed 1010): the one with saved maps.
# It is rebuilt with the simulator that ran the study (commit 64a06b9), since
# later changes to Stage 1 (the resistance path) change the simulated data.
study_code <- tempfile(fileext = ".R")
writeLines(system("git show 64a06b9:R/epiwave-foi-model.R", intern = TRUE), study_code)
study_env <- new.env()
sys.source(study_code, study_env)
d <- study_env$simulate_epiwave_data(seed = 1010L)
stopifnot(isTRUE(all.equal(d$I_true_mat, study$example$truth$I_true_mat)))

thin <- function(m, n = 800) m[round(seq(1, nrow(m), length.out = n)), , drop = FALSE]

replicate1 <- list(
  times = d$times, coords = r4(unname(d$spatial_coords_norm)),
  I_star = r4(d$I_star), I_true = r4(d$I_true_mat),
  epsilon = r4(d$epsilon_true_mat),
  pred_with = r4(study$example$with$pred_mean),
  pred_without = r4(study$example$without$pred_mean),
  cases = d$observed_cases, population = d$true_params$population,
  surveys = d$prev_data[c("survey_indices", "n_tested", "n_positive", "true_prevalence")],
  draws_with = r4(as.data.frame(thin(study$example$with$draws))),
  draws_without = r4(as.data.frame(thin(study$example$without$draws))),
  truth = list(alpha_effective = study$example$truth$true_params$alpha +
                 mean(d$epsilon_true_mat),
               gamma = d$true_params$reporting_rate, sigma2 = d$true_params$gp_sigma^2,
               phi = d$true_params$gp_phi, theta = d$true_params$gp_rho),
  mean_log_I_star = mean(log(pmax(d$I_star, I_STAR_FLOOR))))
rm(d)

study_out <- list(
  meta = study$meta[c("n_reps", "n_sites", "n_times", "n_samples", "warmup", "chains", "date")],
  per_rep = transform(study$per_rep, mean = r4(mean), lwr = r4(lwr), upr = r4(upr),
                      rhat = r4(rhat), ess = r4(ess), true = r4(true)),
  predictive = transform(study$predictive, rmse = r4(rmse), mae = r4(mae),
                         field_coverage = r4(field_coverage)),
  recovery = transform(study$recovery_summary,
                       bias = r4(bias), coverage = r4(coverage),
                       median_rhat = r4(median_rhat), median_ess = r4(median_ess)))

# Reference values for the JavaScript tests
demo <- simulate_epiwave_data()
q_monthly <- transform_convolution_kernel(default_q_daily, MAX_DETECT_DAYS, DAYS_PER_STEP)
set.seed(1)
innov <- matrix(rnorm(3 * 6), 3, 6)
refs <- list(
  q_daily = default_q_daily(0:35),
  q_monthly = q_monthly(0:3),
  ode_site1 = list(I_star = demo$I_star[1, ], x = demo$x_star[1, ], times = demo$times),
  ode_eip10 = local({
    times <- demo$times
    m <- get_fixed_m(times, "s", seasonal_amplitude = 0.6)
    a <- get_fixed_a(times, "s"); g <- get_fixed_g(times, "s")
    adj <- apply_interventions(m, a, g, itn_coverage = outer(0.7, seq(0, 1, length.out = length(times))),
                               itn_susceptibility = 0.8)
    ode <- solve_ross_macdonald_multi_site(adj$m, adj$a, adj$g, times, eip_days = 10)
    list(I_star = as.vector(compute_mechanistic_prediction(adj$m, adj$a, 0.8, ode$z)))
  }),
  ar1 = list(rho = 0.7, innovations = innov, out = ar1(0.7, innov)),
  matern = list(d = c(0, 0.1, 0.5, 1), phi = 0.8,
                k = { r <- sqrt(5) * c(0, 0.1, 0.5, 1) / 0.8; (1 + r + r^2 / 3) * exp(-r) }))

out <- list(replicate1 = replicate1, study = study_out, true_params = TRUE_PARAMS,
            priors = PRIORS, refs = refs)
writeLines(paste0("window.EPIWAVE_DATA = ",
                  jsonlite::toJSON(out, digits = NA, auto_unbox = TRUE, matrix = "rowmajor",
                                   dataframe = "columns"), ";"),
           "site/src/data.js")
cat("wrote site/src/data.js\n")
