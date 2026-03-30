# The Five Parameters: What They Do and Why You Should Care

**EpiWave FOI Model — Stage 2 Parameter Guide**
Ernest Moyo | NM-AIST / Vector Atlas | March 2026

---

## The Setup

Stage 2 of the EpiWave FOI model calibrates the mechanistic prediction I* (from Stage 1) to observed case counts and prevalence surveys. The structural equation is:

    log(I_{s,t}) = alpha + log(I*_{s,t}) + epsilon_{s,t}

    Cases ~ Poisson(gamma * I_{s,t} * N_s)
    Prevalence ~ Binomial(N_tested, x_adjusted)

There are exactly **5 free parameters** estimated via MCMC (HMC in greta/TensorFlow). Everything else — the entomological parameters, the ODE solution, I* itself — is fixed from Stage 1.

This document explains what each parameter controls, what its prior is, what happens if it's wrong, and what questions to ask about it.

---

## 1. alpha — The Intercept

**What it does:** alpha adjusts the overall level of infection incidence relative to the mechanistic prediction I*. It lives on the log scale.

- If alpha = 0: I* is perfectly calibrated. The mechanistic model gets the average incidence right.
- If alpha > 0: true incidence is systematically higher than I* predicts.
- If alpha < 0: true incidence is systematically lower than I* predicts.

**Prior:** Normal(0, 1) — centred at zero (no systematic bias), wide enough to allow substantial correction.

**Why you should care:** alpha tells you how well the mechanistic model (Stage 1) captures average transmission intensity. A large alpha means the entomological parameters from Vector Atlas are systematically off — maybe mosquito abundance is underestimated, or the ODE dynamics don't capture some aspect of transmission. In a real application, alpha absorbs all the systematic errors in the mechanistic model.

**What can go wrong:** alpha is entangled with gamma (reporting rate) through the case count likelihood. For any constant k, shifting alpha up by k and multiplying gamma by exp(-k) gives the same expected case counts. This is the alpha-gamma non-identifiability problem. The prevalence likelihood is supposed to break this, because prevalence data informs alpha independently of gamma.

**Current status:** In our simulation (true alpha = 0), the posterior median is -0.66. This bias is linked to the alpha-gamma ridge — the two parameters are trading off.

---

## 2. gamma (gamma_rr) — The Reporting Rate

**What it does:** gamma is the case ascertainment rate — the fraction of true infections that appear in the routine surveillance system as reported cases.

- gamma = 0.1 means 10% of infections are reported.
- gamma = 0.5 means half of infections are reported.
- gamma = 1.0 means perfect reporting (unrealistic for malaria).

**Prior:** Normal(0.1, 0.05), truncated above 0.001 — centred at 10% reporting, which is typical for malaria in sub-Saharan Africa.

**Why you should care:** gamma directly controls the scale of predicted case counts. If gamma is wrong, your case count predictions are wrong by a multiplicative factor. More importantly, gamma is what connects the biological model (infection incidence) to the operational data (reported cases). Knowing gamma separately from alpha means you can say "the mechanistic model is right about incidence, but only 10% of cases are reported" versus "the mechanistic model overestimates incidence."

**What can go wrong:** gamma is non-identifiable with alpha from case data alone (see alpha above). The prevalence likelihood should resolve this because prevalence surveys measure actual infection in the population, completely independent of whether those infections get reported to the health system. But this requires sufficient prevalence data.

**Current status:** In our simulation (true gamma = 0.1), the posterior median is 0.15. It's trading off with alpha along the ridge.

---

## 3. sigma-squared (sigma2) — GP Marginal Variance

**What it does:** sigma-squared controls how much the true infection incidence departs from the mechanistic prediction I*. It is the variance of the Gaussian Process residuals.

- sigma2 small (e.g. 0.1): the mechanistic model explains most of the spatial variation. The GP makes only minor corrections.
- sigma2 large (e.g. 3.0): the mechanistic model misses a lot. The GP has to do heavy lifting to explain spatial variation that I* doesn't capture.

**Prior:** LogNormal(-0.5, 0.5) — median around 0.6, always positive.

**Why you should care:** sigma2 is the key diagnostic for "how good is my mechanistic model?" If you include the I* offset and sigma2 is small, the entomological data is doing its job — Vector Atlas parameters are explaining spatial variation in malaria burden. If sigma2 is large, something is wrong upstream: maybe the entomological parameters are inaccurate, or there are important spatial drivers (like healthcare access, urbanisation, human mobility) that the Ross-Macdonald model doesn't capture.

The comparison between models is telling: in our simulation, sigma2 = 0.79 WITH the offset vs 3.69 WITHOUT. The mechanistic model reduces the GP's workload by nearly 5x.

**What can go wrong:** sigma2 can absorb bias from alpha. If alpha is wrong (e.g. biased by the ridge), the GP variance inflates to compensate. This is why getting alpha right (via the prevalence likelihood) matters.

**Current status:** True sigma2 = 0.36, estimated at 0.79. Overestimated, likely because alpha bias forces the GP to work harder.

---

## 4. phi — Spatial Lengthscale

**What it does:** phi controls the range of spatial correlation in the GP residuals. It is the lengthscale parameter of the Matern 5/2 kernel.

- phi small (e.g. 0.3): residuals are locally correlated. Nearby sites can have very different departures from I*.
- phi large (e.g. 3.0): residuals are correlated over long distances. If the mechanistic model over-predicts in one location, it likely over-predicts across a wide region.

**Prior:** LogNormal(0.5, 0.5) — median around 1.65, always positive.

**Why you should care:** phi tells you the spatial scale of what the mechanistic model gets wrong. If phi is large, the mechanistic model has broad regional biases — maybe Vector Atlas estimates are systematically off across a whole province. If phi is small, the errors are local — maybe individual sites have unique conditions (irrigation, urbanisation) that the ODE doesn't capture. This has direct implications for how to improve the mechanistic model.

**What can go wrong in our simulation:** Our coordinates are normalised to [0, 1], and the true phi = 3.0. On this scale, Matern 5/2 with phi = 3.0 produces >86% correlation between ALL site pairs. The GP residuals are essentially a shared time-varying signal, not a spatially-varying one. The data cannot distinguish phi = 3.0 from phi = 5.0 or phi = 10.0 — they all look the same. So the prior (median 1.65) pulls the estimate toward 1.5, because that's where the prior mass is and the data offers no resistance.

This is a simulation design issue, not a model problem. With real data on real geographic coordinates (degrees or km), phi would have a meaningful interpretation — e.g. "residuals are correlated over ~200 km."

**Current status:** True phi = 3.0, estimated at 1.53. The data is not informative about phi at this coordinate scale.

---

## 5. theta — AR(1) Temporal Correlation

**What it does:** theta controls how persistent the GP residuals are across time. It is the AR(1) autoregressive coefficient.

- theta = 0: no temporal memory. This month's residual is independent of last month's.
- theta = 0.5: moderate persistence. Residuals decay to half their value each time step.
- theta = 0.9: strong persistence. If the mechanistic model under-predicts this month, it probably under-predicts next month too.

**Prior:** Uniform(0, 1) — no preference within the valid range.

**Why you should care:** theta captures whether the mechanistic model's errors are transient or persistent. High theta means systematic biases that last for months — perhaps the seasonal pattern in mosquito abundance is wrong, or an intervention effect (like ITN distribution) is not accurately modelled. Low theta means the departures are noisy and short-lived.

In combination with phi, theta and phi together describe the full correlation structure of what the mechanistic model gets wrong: phi tells you how far errors extend in space, theta tells you how long they persist in time.

**What can go wrong:** theta is generally the best-behaved parameter. It is directly informed by the temporal autocorrelation in the residuals, and there is no identifiability issue.

**Current status:** True theta = 0.75, estimated at 0.73. This is the best-recovered parameter — confirming that the AR(1) temporal structure is well-identified by the data.

---

## How the Parameters Interact

The five parameters are not independent. Here is how they relate:

    alpha <---> gamma    (non-identifiable from cases alone; prevalence breaks this)
    alpha  --> sigma2    (if alpha is biased, sigma2 inflates to compensate)
    phi   <--> sigma2    (large phi + small sigma2 can look like small phi + large sigma2)
    theta               (largely independent — best identified)

The most important interaction is the alpha-gamma ridge. Everything downstream (sigma2, predictions, maps) depends on getting this right. The dual likelihood (prevalence data) is the mechanism that should resolve it.

---

## Summary Table

| Parameter | Symbol | True | Estimated | Prior | What It Tells You |
|-----------|--------|------|-----------|-------|-------------------|
| Intercept | alpha | 0.0 | -0.66 | N(0, 1) | How well I* captures average incidence |
| Reporting rate | gamma | 0.10 | 0.15 | N(0.1, 0.05)+ | Fraction of infections reported |
| GP variance | sigma2 | 0.36 | 0.79 | LogN(-0.5, 0.5) | How much reality departs from I* |
| Spatial lengthscale | phi | 3.0 | 1.53 | LogN(0.5, 0.5) | Spatial range of mechanistic model errors |
| AR(1) correlation | theta | 0.75 | 0.73 | Uniform(0, 1) | Temporal persistence of errors |

---

## Questions to Bring to Supervisors

1. **Is 30% survey coverage enough to break the alpha-gamma ridge?** Our simulation has 147 prevalence surveys across 490 site-time combinations. The ridge persists. How much prevalence data would be available for the real study site?

2. **Should we adjust the simulation design for phi?** True phi = 3.0 on [0, 1] coordinates is not estimable. Should we use a lower phi, or work in raw coordinate space?

3. **Is the sigma2 overestimate acceptable?** sigma2 = 0.79 vs true 0.36 — is this driven by the alpha bias, or is it a prior issue?

4. **What does theta = 0.75 imply for the study site?** Strong month-to-month persistence in residuals. Does this match expectations about how entomological conditions vary temporally?

---

*This document accompanies the EpiWave FOI Model presentation (March 2026). All numbers are from the baseline simulation-estimation run with 10 sites, 49 time steps, 1000 MCMC samples.*
