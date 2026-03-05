# Issue: Stage 2 GP Approach Replaced with Negative Binomial

**Date:** 2026-03-05
**Author:** Ernest Moyo
**Status:** Resolved
**Related code:** `R/epiwave-foi-model.R` (commit `4bf94df`)
**Experiment:** `R/experiments/gp_stage2_experiment.R`

---

## Summary

Nick's original Stage 2 specification used a Gaussian Process to model residual variation:

```
log(I(t,s)) = alpha + log(I*(t,s)) + epsilon(t,s)
epsilon ~ GP(0, K)
cases ~ Poisson(gamma * I)
```

This was replaced with a Negative Binomial likelihood with a merged log_rate parameter. This document explains why.

## What Was Tried (v1 through v4)

### v1: Joint Space-Time GP
- **Specification:** Single RBF kernel over (lon, lat, time), creating an n_obs x n_obs covariance matrix
- **Problem:** For 10 sites x 49 months = 490 observations, the 490x490 Cholesky decomposition is near-singular
- **Result:** 0% HMC acceptance. Chains completely stuck.

### v2: Separable GP (Space + Time)
- **Specification:** Additive spatial GP + temporal GP to reduce matrix size
- **Problem:** Still ill-conditioned. HMC very slow, poor mixing.
- **Result:** Failed to converge in reasonable time.

### v3: Hierarchical Random Effects
- **Specification:** Dropped GP entirely, used site-level random effects
- **Problem:** alpha and gamma (reporting rate) are non-identifiable — they form a ridge in posterior space. exp(alpha) * gamma is identifiable but the individual parameters are not.
- **Result:** Runs but wastes MCMC budget exploring the non-identifiability ridge. Posteriors misleadingly wide.

### v4: Negative Binomial with merged log_rate (CURRENT)
- **Specification:**
  ```
  log(mu) = log_rate + log(I*)
  cases ~ NegBin(size, size/(size+mu))
  log_rate ~ Normal(-2, 1)
  log_size ~ Normal(3, 1)
  ```
- **Key fixes:**
  1. Merged alpha + log(gamma) into single `log_rate` parameter (resolves non-identifiability)
  2. NegBin handles overdispersion without GP complexity (2 params vs 5+)
  3. `log_size` reparameterisation avoids boundary issues when model is near-Poisson
- **Result:** Clean convergence. ESS 327-489. Rhat < 1.02. Reporting rate recovered: 0.1003 (true = 0.1000).

## Evidence

Run the experiment to reproduce all four versions:
```r
source("R/experiments/gp_stage2_experiment.R")
```

### Simulation-estimation results (v4):

| Metric | WITH offset | WITHOUT offset |
|--------|------------|----------------|
| log_rate mean (sd) | -2.298 (0.018) | -2.283 (0.170) |
| log_size mean (sd) | 5.486 (0.487) | 0.006 (0.214) |
| Interpretation | Near-Poisson (I* explains variance) | Heavy overdispersion (poor fit) |

- 88% RMSE improvement WITH vs WITHOUT offset
- WITH offset: 10x more precise estimate of reporting rate

## Justification for Departure from Original Specification

Nick's GP-based Stage 2 was the right conceptual framework:
- Use I* as an offset (structural prior from vector biology)
- Model residual variation statistically
- Compare WITH vs WITHOUT to demonstrate value of vector data

The NegBin achieves these same goals:
- I* enters as a log-offset (identical role)
- NegBin size parameter captures residual overdispersion (replaces GP variance)
- WITH vs WITHOUT comparison works identically
- Only 2 MCMC parameters instead of 5+ (alpha, gamma, l_space, l_time, sigma2)

The departure is computational, not conceptual. The scientific framework is preserved.

## Future Consideration

If spatial correlation modelling becomes necessary (e.g., for prediction at unobserved locations), options include:
- INLA/SPDE mesh approach (avoids dense GP matrices)
- Nearest-neighbour GP approximation
- Post-hoc spatial smoothing of NegBin residuals

For the simulation-estimation paper (Paper 2), the NegBin is sufficient to demonstrate the framework's value.

## References

- Experiment code: `R/experiments/gp_stage2_experiment.R`
- Production code: `R/epiwave-foi-model.R` lines 726-776
- Commit: `4bf94df` — "Reparameterise Stage 2: log_size ~ Normal(3,1)"
- Documentation: `docs/03_model_evolution/phase4_stage2_iterations.md`
- Architectural decision: `docs/04_key_decisions/architectural_decisions.md` (ADR-004, ADR-005)
