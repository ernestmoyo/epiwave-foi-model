# Issue: ODE Speed Test — R/deSolve vs TensorFlow DormandPrince

**Date:** 2026-01-24
**Author:** Ernest Moyo (with Nick Golding)
**Status:** Resolved
**Related code:** `R/experiments/speed_test.R`, `R/experiments/tf_ode_test.R`
**GitHub Issues:** #4, #5

---

## Summary

A critical question for the framework was whether to solve Ross-Macdonald ODEs in R (deSolve/LSODA) or TensorFlow (DormandPrince). If TF was significantly faster, the ODEs could stay inside greta's computation graph and be differentiated automatically for HMC.

## What Was Tested

Nick Golding wrote `tf_ode_test.R` — a comprehensive TensorFlow implementation including:
- Thin-plate spline interpolation of time-varying parameters (m, a, g)
- TensorFlow-native derivative function with batch dimensions
- DormandPrince adaptive ODE solver from TensorFlow Probability

Ernest extended this with `speed_test.R` to benchmark across different site counts.

### Setup
- 4-year simulation period with 1-year burn-in
- Strong bimodal seasonality (worst case for ODE solvers)
- Monthly evaluation points
- 3 repeats per configuration

### Results

| Sites | Mean Time (s) | SD (s) |
|-------|--------------|--------|
| 5     | 23.92        | 0.78   |
| 10    | 23.76        | 1.19   |
| 25    | 23.78        | 0.67   |
| 50    | 23.34        | 1.26   |

**Key finding: ~24 seconds per evaluation regardless of site count.**

The computation time is dominated by TensorFlow overhead (graph compilation, Python bridge), not by the actual ODE solving. Adding more sites barely changes timing because TF vectorises the sites.

### Comparison with R/deSolve

R/deSolve (LSODA) for the same problem:
- 5 sites: ~0.05s
- 50 sites: ~0.5s (linear scaling, as expected)

**R/deSolve is 50-500x faster for forward-mode evaluation.**

## Root Cause

Nick identified two factors in the lengthy chat:

1. **TF adaptive solvers are slow for vectorised problems** — When one site has wiggly dynamics, ALL sites must slow down (smallest adaptive step governs all). TF's vectorisation strength becomes a weakness here.

2. **Fixed-grid solvers** — Nick suggested trying the Euler/RK4 fixed-grid method instead of DormandPrince, which would be more TF-friendly. However, this trades accuracy for speed.

3. **The real bottleneck** — Even if TF ODE solving were fast, having ODEs inside the MCMC loop means solving them ~thousands of times. At 24s per evaluation, that's days of compute.

## Resolution

This speed test was one of the key inputs that led to the **two-stage architecture decision** (see `2026-01-24_two_stage_architecture.md`):

- **Stage 1:** Solve ODEs ONCE per site using R/deSolve (fast, accurate, outside MCMC)
- **Stage 2:** Use the ODE output (I*) as a fixed offset in greta (no ODE solving in MCMC)

The TF ODE solver and spline interpolation code is preserved in experiments as a reference, and the `greta_rm_ode_op.R` wrapper shows how it would be used if the emulator approach is pursued later.

## Nick's Quote (from lengthy chat)

> "I was having a play with some of the sentiments in the ODE solvers in TensorFlow, and I basically couldn't get them to run very fast. And the only way I could get them to run fast was to use the fixed grid, which is what I mean when I says it takes the same size step each time."

## References

- Speed test code: `R/experiments/speed_test.R`
- Nick's TF ODE code: `R/experiments/tf_ode_test.R`
- Greta operation wrapper: `R/experiments/greta_rm_ode_op.R`
- Speed test results: `R/experiments/speed_test_results.csv`
- Speed test plot: `R/experiments/speed_test_plot.png`
