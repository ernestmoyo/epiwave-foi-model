# Issue: ESS=82 Convergence Failure with Beta(phi) Parameterisation

**Date:** 2026-03-05
**Author:** Ernest Moyo
**Status:** Resolved — reparameterised to log_size
**Related code:** `R/epiwave-foi-model.R` lines 754-760
**Commit:** `4bf94df`

---

## Summary

The initial NegBin Stage 2 model (v4) used a Beta prior on an overdispersion parameter phi, with size = (1-phi)/phi. When the mechanistic offset explained nearly all variation, phi was pushed toward 0 (near-Poisson), creating a boundary issue that caused ESS to collapse to 82.

## The Problem

### Original parameterisation
```r
phi  <- beta(1, 9)            # overdispersion, constrained to (0, 1)
size <- (1 - phi) / phi       # NB size parameter
model(log_rate, phi)
```

### What happened
- WITH mechanistic offset: the model fit was excellent, so overdispersion was minimal
- Minimal overdispersion means size -> infinity, which means phi -> 0
- Beta(1,9) prior near 0 creates a steep gradient that HMC struggles with
- The sampler gets stuck near the boundary: **ESS = 82** (out of 1000 samples)
- Rhat remained acceptable but the chains were barely moving

### Diagnostic output
```
phi: mean=0.0003  sd=0.0002  ESS=82  Rhat=1.04
```

ESS of 82 means only ~8% of samples are effectively independent — the rest are autocorrelated copies. This is far below the recommended minimum of 400.

## The Fix

### New parameterisation
```r
log_size <- normal(3, 1)      # unconstrained real line
size     <- exp(log_size)     # always positive, smooth for HMC
model(log_rate, log_size)
```

### Why this works
1. **No boundary:** log_size lives on (-inf, +inf), so HMC never hits a wall
2. **Smooth geometry:** Normal prior on log scale gives a smooth posterior surface
3. **Near-Poisson is fine:** When size is large (near-Poisson), log_size just gets large — no numerical issues
4. **Prior interpretation:** Normal(3,1) => median size ~20, 95% CI [~1, ~400] — broad but regularising

### Result after fix
```
log_size: mean=5.486  sd=0.487  ESS=489  Rhat=1.004
```

ESS jumped from 82 to **489** — a 6x improvement. Chains mix well, Rhat is excellent.

## Lesson

When a parameter can be pushed to a boundary by the data (phi -> 0 means "no overdispersion needed"), parameterise on an unconstrained scale. The log transform is the standard solution for positive parameters that may need to be very large.

This is a well-known issue in Bayesian computation:
- Stan's documentation recommends the same: use `log(sigma)` not `sigma` when sigma could be near 0
- The "folk theorem of statistical computing" — if your sampler is slow, it's probably a parameterisation problem

## References

- Commit: `4bf94df` — "Reparameterise Stage 2: log_size ~ Normal(3,1) replacing phi ~ Beta(1,9)"
- Production code: `R/epiwave-foi-model.R` lines 754-760
- GP experiment (shows full v1-v4 evolution): `R/experiments/gp_stage2_experiment.R`
- Related non-identifiability fix: `correspondence/issues/2026-03-02_alpha_gamma_nonidentifiability.md`
