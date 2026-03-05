# Issue: Adoption of Two-Stage Architecture (The Brainwave)

**Date:** 2026-01-24
**Author:** Ernest Moyo
**Status:** Resolved — implemented in production code
**Related code:** `R/epiwave-foi-model.R`
**Preceded by:** ODE speed test results (`2026-01-24_ode_speed_test_r_vs_tf.md`)

---

## Summary

The original model design had all parameters (entomological + statistical) estimated jointly inside a single MCMC loop. This required solving the Ross-Macdonald ODEs thousands of times during sampling. Speed tests showed this was computationally infeasible (~24s per TF ODE evaluation).

Nick proposed an alternative: **fix the entomological parameters from Vector Atlas data and solve the ODEs once**, then use the output as a structural offset in a simpler statistical model.

## The Problem

The joint inference approach required:
- m(t,s), a(t,s), g(t,s) as parameters with priors
- ODE solving at every MCMC iteration (inside TensorFlow graph)
- Spline interpolation for time-varying parameters inside the derivative function
- ~10+ parameters to infer simultaneously

At ~24s per ODE evaluation and ~10,000 MCMC iterations needed, this would take **~67 hours** for a single chain on a small problem. Scaling to national level (1000s of pixels) was impossible.

## Nick's Insight (January 2026)

From the lengthy chat, Nick's key statements:

> "Not to do statistical inference on the dynamic part of the model. Instead, use fixed values for entomological parameters (VA estimates or temperature models) and solve dynamics ONCE per pixel to get z_t."

> "We have our maps of abundance, and M is just abundance divided by human populations. That's easy. We can have maps of biting rate quite easily."

> "B, we can generally assume it's spatially constant. We don't have good data, so we can basically ignore it."

The idea: Vector Atlas already produces maps of m, a, g (or will soon). Why estimate what we already know? Use the existing estimates as fixed inputs, solve the ODE once, and focus inference on what we DON'T know — the relationship between mechanistic predictions and observed cases.

## The Two-Stage Solution

### Stage 1: Mechanistic Prediction (no inference)
1. Get fixed m(t,s), a(t,s), g(t,s) from Vector Atlas or temperature models
2. Solve Ross-Macdonald ODE once per site using R/deSolve (~0.01s per site)
3. Compute I*(t,s) = m * a * b * z * N (mechanistic incidence prediction)

### Stage 2: Statistical Offset Model (Bayesian inference)
4. Model: log(mu) = log_rate + log(I*) — I* enters as fixed offset
5. Cases ~ NegBin(size, size/(size+mu))
6. Infer only 2 parameters: log_rate, log_size

### Computational Savings
- **Before:** ~10 params, ODE at every iteration, ~67 hours per chain
- **After:** 2 params, no ODE in MCMC loop, ~12 seconds total

That's a **~20,000x speedup** while preserving the mechanistic framework.

## Why This Works Scientifically

1. **Vector Atlas data exists** — continental-scale maps of mosquito abundance, biting rate, and survival are being produced by the Vector Atlas project. Using them as fixed inputs is not an assumption — it's leveraging existing science.

2. **The offset preserves the biology** — I* carries all the seasonal and spatial structure from the ODE dynamics. The statistical model just needs to scale it (reporting rate) and handle noise (overdispersion).

3. **Separation of concerns** — entomological uncertainty is handled by Vector Atlas. Epidemiological uncertainty is handled by Stage 2. Each team does what they're best at.

4. **Nick's analogy from the chat:** The mechanistic prediction I* is like the vectorial capacity — it collapses all the vector biology into one number per location per time. The statistical model then uses this as "beta" in what becomes essentially an SIS-type observation model.

## Implementation

The two-stage architecture is fully implemented in `R/epiwave-foi-model.R`:
- Stage 1 functions: `get_fixed_m()`, `get_fixed_a()`, `get_fixed_g()`, `apply_interventions()`, `solve_ross_macdonald_multi_site()`, `compute_mechanistic_prediction()`
- Stage 2 function: `fit_epiwave_with_offset()`
- Validation: `simulate_and_estimate()`

## References

- Production code: `R/epiwave-foi-model.R`
- Speed test that motivated this: `R/experiments/speed_test.R`
- Nick's TF ODE work (preserved): `R/experiments/tf_ode_test.R`
- Meeting transcript: `meetings/PhD_Work/Lengthy Chat with Nick about Code.txt`
- Architecture decisions: `docs/04_key_decisions/architectural_decisions.md`
