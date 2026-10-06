# ==============================================================================
# Publication-grade simulation-estimation harness
# ==============================================================================
# Runs many simulated replicates, fits both models (WITH I* offset vs I*=0), and
# scores parameter recovery, credible-interval coverage, MCMC convergence
# (Rhat / ESS), and predictive accuracy of the latent incidence surface.
#
# Depends on R/epiwave-foi-model.R (sourced) and a working greta install
# (source R/greta_setup.R first). Does NOT load ggplot2, so it runs headless.
#
# Usage:
#   Sys.setenv(RETICULATE_PYTHON = "C:/Users/ernes/AppData/Local/r-miniconda/envs/r-greta/python.exe")
#   source("R/greta_setup.R"); source("R/epiwave-foi-model.R")
#   source("R/sim_estimation_harness.R")
#   study <- run_sim_estimation_study(n_reps = 50, chains = 4,
#                                     n_samples = 2000, warmup = 2000)
# ==============================================================================

PARAMS <- c("alpha", "gamma_rr", "sigma2", "phi", "theta")

#' Per-parameter posterior summary + recovery + convergence for one model fit.
#' @keywords internal
.score_parameters <- function(draws, truth_vec, center_offset) {
  mat <- as.matrix(draws)
  rh  <- tryCatch(coda::gelman.diag(draws, multivariate = FALSE)$psrf[, 1],
                  error = function(e) setNames(rep(NA_real_, length(PARAMS)), PARAMS))
  ess <- tryCatch(coda::effectiveSize(draws),
                  error = function(e) setNames(rep(NA_real_, length(PARAMS)), PARAMS))
  truth <- truth_vec
  truth["alpha"] <- truth["alpha"] + center_offset   # effective intercept if centred
  do.call(rbind, lapply(PARAMS, function(p) {
    x <- mat[, p]
    lwr <- unname(quantile(x, 0.025)); upr <- unname(quantile(x, 0.975))
    data.frame(param = p, true = unname(truth[p]),
               mean = mean(x), lwr = lwr, upr = upr,
               covered = truth[p] >= lwr & truth[p] <= upr,
               rhat = unname(rh[p]), ess = unname(ess[p]))
  }))
}

#' Posterior incidence prediction + accuracy for one model fit.
#' @keywords internal
.score_prediction <- function(model, draws, I_true_mat, held_out = NULL) {
  I_node <- attr(model, "I_latent")
  ns <- nrow(I_true_mat); nt <- ncol(I_true_mat)
  cc <- as.matrix(greta::calculate(I_node, values = draws))   # [n_draws x (ns*nt)], column-major
  pred_mean <- matrix(colMeans(cc), ns, nt)
  pred_lwr  <- matrix(apply(cc, 2, quantile, 0.025), ns, nt)
  pred_upr  <- matrix(apply(cc, 2, quantile, 0.975), ns, nt)
  score <- compute_performance_metrics(pred_mean, I_true_mat, pred_lwr, pred_upr)
  # map error at sites the fit never saw (out-of-sample), when there are any
  held <- if (length(held_out)) compute_performance_metrics(
    pred_mean[held_out, ], I_true_mat[held_out, ], pred_lwr[held_out, ], pred_upr[held_out, ]) else
    list(rmse = NA, coverage = NA)
  list(rmse = score$rmse, mae = score$mae, field_coverage = score$coverage,
       rmse_heldout = held$rmse, coverage_heldout = held$coverage,
       pred_mean = pred_mean)
}

#' Run the full multi-replicate simulation-estimation study.
#'
#' @param n_reps Number of simulated replicates
#' @param n_sites,n_times Grid dimensions
#' @param n_samples,warmup,chains MCMC settings (>=2 chains for Rhat)
#' @param center_alpha Centre the GP field (improves alpha mixing)
#' @param use_sparse_gp,n_inducing Sparse GP options (passed to fit)
#' @param true_params Optional true-parameter list (NULL = defaults)
#' @param holdout_fraction Share of sites whose data are withheld from the fit
#'   and used only to score the map (0 = fit and score everything)
#' @param site_variation Give sites different ITN coverage and susceptibility,
#'   so I* varies in space (see simulate_epiwave_data)
#' @param base_seed Seed of replicate 1 (replicate r uses base_seed + 10*r)
#' @param out_path Where to save the study object (.RData)
#' @return A `study` list (also saved to out_path)
#' @export
run_sim_estimation_study <- function(n_reps = 50,
                                     n_sites = 10, n_times = 48,
                                     n_samples = 2000, warmup = 2000, chains = 4,
                                     center_alpha = TRUE,
                                     use_sparse_gp = FALSE, n_inducing = 40,
                                     true_params = TRUE_PARAMS,
                                     holdout_fraction = 0, site_variation = FALSE,
                                     base_seed = 1000,
                                     out_path = "outputs/sim_estimation_results.RData") {

  stopifnot(chains >= 2)
  fit_one <- function(d, use_mech, observed_sites) {
    inducing <- if (use_sparse_gp && n_inducing < n_sites)
      d$spatial_coords_norm[round(seq(1, n_sites, length.out = n_inducing)), , drop = FALSE] else NULL
    m <- fit_epiwave_gp(observed_cases = d$observed_cases, I_star = d$I_star,
                        N_pop = d$pop_matrix, spatial_coords = d$spatial_coords_norm,
                        prev_data = d$prev_data, prev_conv_matrix = d$conv_matrix,
                        use_mechanistic = use_mech, inducing = inducing,
                        center_alpha = center_alpha, observed_sites = observed_sites)
    draws <- greta::mcmc(m, n_samples = n_samples, warmup = warmup,
                         chains = chains, verbose = FALSE)
    list(model = m, draws = draws)
  }

  per_rep <- list(); pred_rows <- list(); example <- NULL

  for (r in seq_len(n_reps)) {
    message(sprintf("\n===== replicate %d / %d =====", r, n_reps))
    d <- simulate_epiwave_data(n_sites = n_sites, n_times = n_times,
                               true_params = true_params, seed = base_seed + 10L * r,
                               site_variation = site_variation)
    # held-out sites, drawn with their own seed so the data are unchanged
    held_out <- if (holdout_fraction > 0) {
      set.seed(base_seed + 10L * r + 4L)
      sort(sample(n_sites, max(1, round(holdout_fraction * n_sites))))
    } else integer(0)
    observed_sites <- if (length(held_out)) setdiff(seq_len(n_sites), held_out) else NULL
    tp <- d$true_params
    truth_vec <- c(alpha = tp$alpha, gamma_rr = tp$reporting_rate,
                   sigma2 = tp$gp_sigma^2, phi = tp$gp_phi, theta = tp$gp_rho)
    center_offset <- if (center_alpha) mean(d$epsilon_true_mat) else 0

    for (mod in c("with", "without")) {
      fit <- tryCatch(fit_one(d, use_mech = (mod == "with"), observed_sites = observed_sites),
                      error = function(e) { message("  fit error (", mod, "): ", e$message); NULL })
      if (is.null(fit) || is.null(fit$draws)) next
      sp <- .score_parameters(fit$draws, truth_vec, center_offset)
      sp$rep <- r; sp$model <- mod
      per_rep[[length(per_rep) + 1]] <- sp

      pr <- .score_prediction(fit$model, fit$draws, d$I_true_mat, held_out)
      pred_rows[[length(pred_rows) + 1]] <- data.frame(
        rep = r, model = mod, rmse = pr$rmse, mae = pr$mae,
        field_coverage = pr$field_coverage,
        rmse_heldout = pr$rmse_heldout, coverage_heldout = pr$coverage_heldout)

      if (r == 1) {
        # Store plain matrices (not the live greta_mcmc_list) so the saved
        # RData is portable and loads without greta attached.
        example[[mod]] <- list(draws = as.matrix(fit$draws), pred_mean = pr$pred_mean)
      }
    }
    if (r == 1) {
      example$truth <- list(I_true_mat = d$I_true_mat, prev_data = d$prev_data,
                            epsilon_true_mat = d$epsilon_true_mat, true_params = tp)
    }
  }

  per_rep <- do.call(rbind, per_rep)
  predictive <- do.call(rbind, pred_rows)

  # ---- aggregate summaries ----
  recovery_summary <- do.call(rbind, lapply(split(per_rep, list(per_rep$model, per_rep$param)),
    function(s) data.frame(
      model = s$model[1], param = s$param[1], true = s$true[1],
      bias = mean(s$mean - s$true), rmse_pointest = sqrt(mean((s$mean - s$true)^2)),
      mean_ci_width = mean(s$upr - s$lwr), coverage = mean(s$covered),
      median_rhat = median(s$rhat, na.rm = TRUE), median_ess = median(s$ess, na.rm = TRUE),
      n_rep = nrow(s))))
  rownames(recovery_summary) <- NULL

  predictive_summary <- do.call(rbind, lapply(split(predictive, predictive$model),
    function(s) data.frame(model = s$model[1],
      mean_rmse = mean(s$rmse), mean_mae = mean(s$mae),
      mean_field_coverage = mean(s$field_coverage),
      mean_rmse_heldout = mean(s$rmse_heldout), mean_coverage_heldout = mean(s$coverage_heldout),
      n_rep = nrow(s))))
  rownames(predictive_summary) <- NULL

  study <- list(
    meta = list(n_reps = n_reps, n_sites = n_sites, n_times = n_times,
                n_samples = n_samples, warmup = warmup, chains = chains,
                center_alpha = center_alpha, use_sparse_gp = use_sparse_gp,
                holdout_fraction = holdout_fraction, site_variation = site_variation,
                true_params = true_params,
                priors = PRIORS,
                git_sha = tryCatch(system("git describe --always --dirty", intern = TRUE),
                                   error = function(e) NA_character_),
                date = as.character(Sys.Date())),
    per_rep = per_rep, predictive = predictive,
    recovery_summary = recovery_summary, predictive_summary = predictive_summary,
    example = example)

  if (!is.null(out_path)) {
    dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
    save(study, file = out_path)
    message("\nSaved study to ", out_path)
  }
  study
}
