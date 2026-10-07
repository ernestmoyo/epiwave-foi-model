# Checks for the pieces borrowed from the Vector Atlas models (2026-10-06).
#   Rscript tests/stage1_checks.R      (from the project root; no greta needed)

suppressPackageStartupMessages(source("R/epiwave-foi-model.R"))

failed <- 0
check <- function(name, ok) {
  if (!isTRUE(ok)) failed <<- failed + 1
  cat(if (isTRUE(ok)) "PASS " else "FAIL ", name, "\n")
}

base <- list(m = matrix(2, 1, 1), a = matrix(0.3, 1, 1), g = matrix(0.1, 1, 1))
itn <- function(...) do.call(apply_interventions, c(base, list(itn_coverage = 0.6, ...)))

# Resistance: full susceptibility keeps the whole ITN effect ...
full <- itn(itn_susceptibility = 1)
check("full susceptibility: m reduced by coverage x kill rate",
      isTRUE(all.equal(full$m[1], 2 * (1 - 0.6 * 0.5), tolerance = 1e-12)))
# ... and none still keeps 54% of it (effective coverage 0.6 x 0.54)
none <- itn(itn_susceptibility = 0)
check("zero susceptibility keeps 54% of the ITN effect",
      isTRUE(all.equal(c(none$m, none$a, none$g),
                       c(2 * (1 - 0.324 * 0.5), 0.3 * (1 - 0.324 * 0.3), 0.1 * (1 + 0.324 * 0.3)),
                       tolerance = 1e-12)))
check("itn_effect_retained(0.8) = 0.908", isTRUE(all.equal(itn_effect_retained(0.8), 0.908)))

# IRS is not affected by ITN (pyrethroid) susceptibility
irs_a <- apply_interventions(base$m, base$a, base$g, irs_coverage = 0.5, itn_susceptibility = 0)
irs_b <- apply_interventions(base$m, base$a, base$g, irs_coverage = 0.5, itn_susceptibility = 1)
check("IRS ignores ITN susceptibility", identical(irs_a, irs_b))
irs_r <- apply_interventions(base$m, base$a, base$g, irs_coverage = 0.5, irs_susceptibility = 0)
check("IRS with zero susceptibility has no effect", identical(irs_r$g, base$g))

# m from the interim biting-rate maps
ma <- matrix(c(0.6, 1.2), 1)
check("m_from_biting_rate inverts m*a",
      isTRUE(all.equal(m_from_biting_rate(ma, matrix(0.3, 1, 2)), matrix(c(2, 4), 1))))

# EIP: NULL leaves the model exactly as before
times <- seq(0, 48 * 30, by = 30)
m1 <- get_fixed_m(times, "s1", seasonal_amplitude = 0.6)
a1 <- get_fixed_a(times, "s1"); g1 <- get_fixed_g(times, "s1")
no_eip <- solve_ross_macdonald_multi_site(m1, a1, g1, times)
null_eip <- solve_ross_macdonald_multi_site(m1, a1, g1, times, eip_days = NULL)
check("eip_days = NULL is identical to the 2-state model", identical(no_eip, null_eip))

# EIP: survival through k exposed stages is (1 + g n / k)^-k, so R0 scales by it
check("eip survival factor, k = 4", isTRUE(all.equal(eip_survival(g = 0.1, eip_days = 10, stages = 4), (1 + 0.25)^-4)))
check("eip survival factor tends to exp(-g n)", abs(eip_survival(0.1, 10, 400) - exp(-1)) < 1e-3)
with_eip <- solve_ross_macdonald_multi_site(m1, a1, g1, times, eip_days = 10)
I_no <- compute_mechanistic_prediction(m1, a1, 0.8, no_eip$z)
I_eip <- compute_mechanistic_prediction(m1, a1, 0.8, with_eip$z)
check("a 10-day EIP lowers I*", mean(I_eip) < mean(I_no))
cat(sprintf("     mean I* without EIP %.3f, with 10-day EIP %.3f (ratio %.2f)\n",
            mean(I_no), mean(I_eip), mean(I_no) / mean(I_eip)))
check("eip_gething: 111 / (T - 16), Inf at or below 16 C",
      isTRUE(all.equal(eip_gething(c(15, 16, 27)), c(Inf, Inf, 111 / 11))))
cold <- solve_ross_macdonald_multi_site(m1, a1, g1, times, eip_days = Inf)
check("no sporogony below 16 C: mosquitoes never become infectious", max(cold$z) < 1e-3 * 1.0001)

# Spatial signal: site_variation gives I* that differs between sites
flat <- simulate_epiwave_data(n_sites = 6, n_times = 12, seed = 3L)
vary <- simulate_epiwave_data(n_sites = 6, n_times = 12, seed = 3L, site_variation = TRUE)
check("default simulation: I* identical across sites", max(apply(flat$I_star, 2, sd)) == 0)
check("site_variation: I* differs across sites", min(apply(vary$I_star[, -1], 2, sd)) > 0)

# Equilibrium start: the ODE started at rm_equilibrium() stays there when the
# parameters are held constant (with and without the EIP stages)
yr <- seq(0, 360, by = 30)
const <- function(v) matrix(v, 1, length(yr))
for (eip in list(NULL, 10)) {
  eq <- rm_equilibrium(m = 0.3, a = 0.3, g = 0.1, b = 0.5, c = 0.5, r = 1/180,
                       eip_days = eip, stages = 4)
  run <- solve_ross_macdonald_multi_site(const(0.3), const(0.3), const(0.1), yr,
                                         b = 0.5, c = 0.5, r = 1/180,
                                         eip_days = eip, start = eq)
  check(sprintf("rm_equilibrium is a fixed point (EIP %s)", if (is.null(eip)) "off" else eip),
        max(abs(run$x / eq$x - 1), abs(run$z / eq$z - 1)) < 1e-6)
}

# Literature ranges (docs/2026-10-07_parameter_ranges)
tp <- TRUE_PARAMS
eq0 <- rm_equilibrium(tp$baseline_m, tp$baseline_a, tp$baseline_g, tp$b, tp$c, tp$r, eip_days = tp$eip_days)
check(sprintf("prevalence before nets is 0.48 (got %.3f)", eq0$x), abs(eq0$x - 0.48) < 0.01)
sim <- simulate_epiwave_data(seed = 1010L)
cat(sprintf("     I*: mean %.4f per day (%.2f per year), range %.1e to %.1e; ODE prevalence %.2f -> %.2f
",
            mean(sim$I_star), 365 * mean(sim$I_star), min(sim$I_star), max(sim$I_star),
            mean(sim$x_star[, 1]), mean(sim$x_star[, ncol(sim$x_star)])))
check("simulated I* is a realistic rate (0.001 to 0.01 per day on average)",
      mean(sim$I_star) > 0.001 && mean(sim$I_star) < 0.01)
check("I* never falls to the floor", min(sim$I_star) > 10 * I_STAR_FLOOR)
exp_cases <- sim$true_params$reporting_rate * sim$I_true_mat * sim$pop_matrix * DAYS_PER_STEP
cat(sprintf("     expected cases per site-month: median %.0f
", median(exp_cases)))
check("expected cases per site-month are informative (median > 20)", median(exp_cases) > 20)

if (failed) { cat("\n", failed, "check(s) failed\n"); quit(status = 1) }
cat("\nAll Stage 1 checks passed.\n")
