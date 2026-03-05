# GitHub Issue #5: TensorFlow ODE Solver Integration

**GitHub:** https://github.com/ernestmoyo/epiwave-foi-model/issues/5
**Date opened:** ~2026-01-24
**Status:** Resolved (code written by Nick, superseded by two-stage approach)
**Related commit:** `7f35603`

---

## Original Request

Integrate TensorFlow Probability's ODE solvers into the greta model so that Ross-Macdonald dynamics can be solved inside the MCMC computation graph, enabling automatic differentiation for HMC.

## What Was Built

Nick Golding wrote `R/experiments/tf_ode_test.R` — a comprehensive implementation including:

1. **Spline interpolation in TensorFlow** — Thin-plate RBF splines to evaluate time-varying parameters at arbitrary time points within the ODE solver
2. **Batched TF derivative function** — Ross-Macdonald equations operating on tensors with batch dimensions for greta compatibility
3. **DormandPrince solver wrapper** — `integrate_RMd()` function calling TFP's adaptive solver
4. **Multi-site vectorisation** — Simultaneous solving for N sites

Ernest then created:
- `R/experiments/speed_test.R` — Benchmarking the TF solver across site counts
- `R/experiments/greta_rm_ode_op.R` — Wrapping the TF solver as a greta operation

## Outcome

The code works correctly — TF ODE solutions match R/deSolve solutions. However, performance was ~24s per evaluation regardless of site count, making it impractical for MCMC (which requires thousands of evaluations).

Nick noted at the bottom of `tf_ode_test.R`:
```r
# this needs debugging, and probably an input signature to be efficient!
```

The `tf_function` compilation step was never fully optimised, which may account for some of the overhead.

## Resolution

Superseded by the two-stage architecture. ODE solving moved to R/deSolve (Stage 1) outside of MCMC. The TF code is preserved in `R/experiments/` as:
- A reference for how to solve ODEs in TensorFlow with spline interpolation
- A starting point for the neural network emulator approach (Paper 2 roadmap)

## Related Files
- Nick's TF code: `R/experiments/tf_ode_test.R`
- Speed test: `R/experiments/speed_test.R`
- Greta wrapper: `R/experiments/greta_rm_ode_op.R`
- Correspondence: `correspondence/issues/2026-01-24_ode_speed_test_r_vs_tf.md`
