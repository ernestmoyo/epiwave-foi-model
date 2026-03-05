# GitHub Issue #4: ODE Solver Speed Benchmarks

**GitHub:** https://github.com/ernestmoyo/epiwave-foi-model/issues/4
**Date opened:** ~2026-01-24
**Status:** Resolved (partially completed, then superseded by two-stage architecture)
**Related commit:** `7f35603`

---

## Original Request

Benchmark different ODE solvers for the Ross-Macdonald model to find the most computationally efficient method for:
- Time-varying parameters (bimodal seasonality)
- Multiple sites (N = 10, 50, 100, 500)
- 4 years of simulation
- Support for automatic differentiation (for Bayesian inference)

## Solvers Tested

### R (deSolve)
- `lsoda` — Adaptive (FORTRAN, stiff/non-stiff switching)
- `rk4` — Fixed-step Runge-Kutta 4th order
- `euler` — Fixed-step Euler
- `ode45` — Dormand-Prince adaptive

### TensorFlow Probability
- `DormandPrince` — Adaptive RK4(5)
- `BDF` — Backward Differentiation Formula (stiff systems)
- Custom fixed-step RK4 and Euler

## Key Finding

TF DormandPrince: ~24s per evaluation, constant regardless of site count (5-50 sites).
R deSolve LSODA: ~0.01s per site, linear scaling.

R was 50-500x faster for forward evaluation.

## Task Checklist (from issue)
- [x] Set up Ross-MacDonald ODEs with time-varying parameters
- [x] Write R baseline tests with deSolve
- [x] Write TensorFlow code for ODE solving
- [x] Speed-test forward-mode solution with different solvers
- [ ] Speed-test differentiation with different solvers (superseded)
- [ ] Wrap optimal solver as greta operation (superseded)

## Resolution

The speed test results directly motivated the **two-stage architecture** decision — solving ODEs in R outside of MCMC rather than inside TensorFlow. The remaining checklist items (differentiation speed test, greta operation wrapping) became unnecessary once we adopted the two-stage approach.

The greta operation wrapper was still created (`R/experiments/greta_rm_ode_op.R`) as a reference for the emulator approach discussed for Paper 2.

## Related Files
- `R/ode_speed_test/` (deleted from main, was committed in wrong branch)
- `R/experiments/speed_test.R` (current location)
- `R/experiments/tf_ode_test.R` (Nick's code)
- Correspondence: `correspondence/issues/2026-01-24_ode_speed_test_r_vs_tf.md`
