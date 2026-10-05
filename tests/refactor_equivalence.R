# Proves a refactor of R/epiwave-foi-model.R changed no numbers.
#
#   Rscript tests/refactor_equivalence.R freeze <model_file>   # save the reference
#   Rscript tests/refactor_equivalence.R check                  # compare the current code
#
# Fingerprints: the simulated data (demo seeds, harness replicate 1, a small
# no-intervention case) and, for five Stage 2 configurations, the size of the
# free state and the joint log-density at three pinned points.

Sys.setenv(RETICULATE_PYTHON = "C:/Users/ernes/AppData/Local/r-miniconda/envs/r-greta/python.exe")
args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]
model_file <- if (mode == "freeze") args[2] else "R/epiwave-foi-model.R"
ref_path <- "cache/refactor_reference_2026-10-05.rds"

suppressPackageStartupMessages({
  source("R/greta_setup.R")
  library(greta)
  source(model_file)
})

# Older code kept the inputs in [n_times x n_sites]; the refactor keeps
# everything in [n_sites x n_times]. Compare in one orientation.
as_sites_by_times <- function(d) {
  n_sites <- nrow(d$spatial_coords_norm)
  flip <- function(x) if (is.matrix(x) && ncol(x) == n_sites && nrow(x) != n_sites) t(x) else x
  for (nm in c("observed_cases", "I_star", "x_star", "pop_matrix")) d[[nm]] <- flip(d[[nm]])
  d
}

simulated <- list(
  demo  = simulate_epiwave_data(),
  rep1  = simulate_epiwave_data(seed = 1010L),
  noint = simulate_epiwave_data(n_sites = 4, n_times = 12,
                                include_interventions = FALSE, seed = 7L))
simulated <- lapply(simulated, as_sites_by_times)

log_density_at <- function(model, points) {
  lp <- model$dag$generate_log_prob_function(which = "unadjusted")
  vapply(points, function(p) {
    fs <- tensorflow::tf$constant(matrix(p, nrow = 1), dtype = tensorflow::tf$float64)
    as.numeric(lp(fs))
  }, numeric(1))
}

d <- simulate_epiwave_data(n_sites = 6, n_times = 12, seed = 5L)
configs <- list(
  with_offset        = list(use_mechanistic = TRUE,  center_alpha = FALSE),
  without_offset     = list(use_mechanistic = FALSE, center_alpha = FALSE),
  with_centred       = list(use_mechanistic = TRUE,  center_alpha = TRUE),
  without_centred    = list(use_mechanistic = FALSE, center_alpha = TRUE),
  with_inducing      = list(use_mechanistic = TRUE,  center_alpha = TRUE,
                            inducing = d$spatial_coords_norm[c(1, 3, 5), ]))

stage2 <- lapply(configs, function(cfg) {
  model <- do.call(fit_epiwave_gp, c(list(
    observed_cases = d$observed_cases, I_star = d$I_star, N_pop = d$pop_matrix,
    spatial_coords = d$spatial_coords_norm, prev_data = d$prev_data,
    prev_conv_matrix = d$conv_matrix), cfg))
  n_free <- sum(lengths(model$dag$example_parameters(free = TRUE)))
  points <- lapply(c(-0.4, 0, 0.3), function(v) {
    set.seed(round(100 * v) + 50)
    v + rnorm(n_free, 0, 0.1)
  })
  list(n_free = n_free, log_density = log_density_at(model, points))
})

fingerprint <- list(simulated = simulated, stage2 = stage2)

if (mode == "freeze") {
  dir.create(dirname(ref_path), showWarnings = FALSE)
  saveRDS(fingerprint, ref_path)
  cat("Reference frozen from", model_file, "\n")
  print(sapply(stage2, function(s) c(n_free = s$n_free, s$log_density)))
} else {
  ref <- readRDS(ref_path)
  ok_sim <- mapply(function(a, b) isTRUE(all.equal(a[names(b)], b, tolerance = 0)),
                   fingerprint$simulated, ref$simulated)
  ok_s2 <- mapply(function(a, b) a$n_free == b$n_free &&
                    isTRUE(all.equal(a$log_density, b$log_density, tolerance = 1e-10)),
                  fingerprint$stage2, ref$stage2)
  print(ok_sim); print(ok_s2)
  if (!all(ok_sim, ok_s2)) {
    for (k in names(ok_sim)[!ok_sim]) print(all.equal(fingerprint$simulated[[k]][names(ref$simulated[[k]])], ref$simulated[[k]]))
    stop("Refactor changed the numbers.")
  }
  cat("PASS: simulated data identical; Stage 2 log-density identical in all configs.\n")
}
