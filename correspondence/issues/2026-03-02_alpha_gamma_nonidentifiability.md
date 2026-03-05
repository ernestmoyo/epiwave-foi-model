# Issue: Alpha/Gamma Non-Identifiability in Stage 2

**Date:** 2026-03-02
**Author:** Ernest Moyo
**Status:** Resolved — merged into single `log_rate` parameter
**Related code:** `R/epiwave-foi-model.R` lines 749-752
**Experiment:** `R/experiments/gp_stage2_experiment.R` (v3)

---

## Summary

When Stage 2 had separate intercept (alpha) and reporting rate (gamma) parameters, MCMC chains explored a non-identifiability ridge instead of converging to a point estimate. Only the product exp(alpha) * gamma is identifiable from the data.

## The Model (before fix)

```
log(mu) = alpha + log(I*) + epsilon
cases ~ Poisson(gamma * mu)
```

Expanding: E[cases] = gamma * exp(alpha) * I* * exp(epsilon)

The data only informs gamma * exp(alpha) — the product. You can increase alpha and decrease gamma by the same factor and get identical predictions. This is a classic non-identifiability.

## Symptoms

- MCMC trace plots showed alpha and gamma anti-correlated (banana-shaped posterior)
- Wide marginal posteriors for both parameters individually
- Product exp(alpha) * gamma was well-estimated, but individual params were not
- Wasted MCMC budget exploring the ridge instead of sampling efficiently

## The Fix

Merge alpha and gamma into a single parameter:

```
log_rate = log(gamma) + alpha
```

So the model becomes:
```
log(mu) = log_rate + log(I*)
cases ~ NegBin(size, size/(size + mu))
```

Now `log_rate` is directly identifiable:
- exp(log_rate) = reporting_rate (the scaling between I* and observed cases)
- Prior: log_rate ~ Normal(-2, 1) centres on exp(-2) ~ 13% reporting

## Validation

Simulation with true reporting rate = 0.1 (log(0.1) = -2.302):
- Recovered: log_rate = -2.298 (sd = 0.018)
- exp(-2.298) = 0.1003 — essentially perfect recovery

## Lesson

When two parameters enter a model only as a product, they are fundamentally non-identifiable. No amount of data or better priors can fix this — you must reparameterise to estimate only the identifiable combination.

## References

- Experiment showing the ridge: `R/experiments/gp_stage2_experiment.R` (v3 section)
- Production fix: `R/epiwave-foi-model.R` lines 749-752
- Related: `docs/04_key_decisions/architectural_decisions.md` (ADR-005)
