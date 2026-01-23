
# =============================================================================
# SPEED TEST: TensorFlow ODE Solver - Varying Number of Sites
# Author: Ernest Moyo
# Building on Nick Golding's integrate_RMd function (PR #5)
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  TensorFlow ODE Solver Speed Test - Ross-MacDonald Model\n
")
cat("================================================================\n\n")

# Ensure Nick's script has been sourced first
if (!exists("integrate_RMd")) {
  cat("Sourcing Nick's tf_ode_test.R...\n")
  source("R/experiments/tf_ode_test.R")
}

# Fixed parameters
params <- list(b = 0.8, c = 0.8, r = 1/7)
days_evaluate <- seq(0, 4*365, by = 30)
burnin_days <- 365

# Function to run speed test for N sites
run_tf_speed_test <- function(n_sites, n_repeats = 3, verbose = TRUE) {
  
  if (verbose) cat(sprintf("  Testing N = %d sites...\n", n_sites))
  
  # Setup data
  previous_days_observed <- -rev(seq(0, burnin_days, by = 30)[-1])
  days_observed <- c(previous_days_observed, days_evaluate)
  
  # Create multi-site data
  set.seed(42)
  m_star_observed <- t(replicate(n_sites, m(days_observed) * exp(rnorm(1, 0, 0.1))))
  a_star_observed <- t(replicate(n_sites, a(days_observed)))
  g_star_observed <- t(replicate(n_sites, g(days_observed)))
  
  # Knots

  n_knots <- length(days_observed)
  epsilon <- 1e-3
  limits <- c(-burnin_days - epsilon, max(days_evaluate) + epsilon)
  knots <- define_knots(n_knots, limits)
  
  # Coefficient matrices
  m_star_coefs <- get_spline_coefs(log(m_star_observed), days_observed, knots)
  a_star_coefs <- get_spline_coefs(qlogis(a_star_observed), days_observed, knots)
  g_star_coefs <- get_spline_coefs(log(g_star_observed), days_observed, knots)
  
  # Timing
  times <- numeric(n_repeats)
  
  for (i in seq_len(n_repeats)) {
    start_time <- Sys.time()
    
    integral <- integrate_RMd(
      burnin = burnin_days,
      times = days_evaluate,
      x_0 = 0.1, z_0 = 0.1,
      m_int = as_tensor(array(0, dim = c(1, 1, 1))),
      m_slope = as_tensor(array(1, dim = c(1, 1, 1))),
      a_int = as_tensor(array(0, dim = c(1, 1, 1))),
      a_slope = as_tensor(array(1, dim = c(1, 1, 1))),
      g_int = as_tensor(array(0, dim = c(1, 1, 1))),
      g_slope = as_tensor(array(1, dim = c(1, 1, 1))),
      b = params$b, c = params$c, r = params$r,
      knots = knots,
      m_star_coefs = m_star_coefs,
      a_star_coefs = a_star_coefs,
      g_star_coefs = g_star_coefs
    )
    
    # Force evaluation
    result <- as.array(integral)
    
    end_time <- Sys.time()
    times[i] <- as.numeric(difftime(end_time, start_time, units = "secs"))
  }
  
  if (verbose) cat(sprintf("    Mean: %.3f s (SD: %.3f)\n", mean(times), sd(times)))
  
  data.frame(
    n_sites = n_sites,
    mean_time = mean(times),
    sd_time = sd(times),
    min_time = min(times),
    max_time = max(times)
  )
}

# Run tests
cat("Running speed tests...\n\n")

n_sites_list <- c(5, 10, 25, 50)
results_list <- list()

for (i in seq_along(n_sites_list)) {
  results_list[[i]] <- run_tf_speed_test(n_sites_list[i], n_repeats = 3)
}

results <- do.call(rbind, results_list)

cat("\n================================================================\n")
cat("  RESULTS SUMMARY\n")
cat("================================================================\n\n")
print(results)

# Save results
write.csv(results, "R/experiments/speed_test_results.csv", row.names = FALSE)
cat("\nResults saved to: R/experiments/speed_test_results.csv\n")

# Plot
png("R/experiments/speed_test_plot.png", width = 800, height = 500)
par(mar = c(5, 5, 4, 2))
plot(results$n_sites, results$mean_time,
     type = "b", pch = 19, col = "steelblue", lwd = 2,
     xlab = "Number of Sites (N)",
     ylab = "Computation Time (seconds)",
     main = "TensorFlow ODE Solver Speed Test\nRoss-MacDonald with DormandPrince",
     cex.lab = 1.2, cex.main = 1.3)
arrows(results$n_sites, results$mean_time - results$sd_time,
       results$n_sites, results$mean_time + results$sd_time,
       angle = 90, code = 3, length = 0.1, col = "steelblue")
grid()
dev.off()
cat("Plot saved to: R/experiments/speed_test_plot.png\n")

cat("\n================================================================\n")
cat("  SPEED TEST COMPLETE\n")
cat("================================================================\n")

