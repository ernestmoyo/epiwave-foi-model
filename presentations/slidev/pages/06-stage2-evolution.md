---
layout: default
---

# Stage 2 Evolution: Learning from Each Iteration

Six iterations — each failure taught us something, converging on the epiwave.mapping-aligned design

<div class="grid grid-cols-3 gap-2 mt-4">

<div class="p-2 rounded-lg bg-red-50 text-xs border border-red-300">

### <span class="text-red-600">v1–v3</span>
**Early GP attempts**

v1: Joint 490x490 GP — Cholesky singular

v2: Additive separable — ill-conditioned

v3: Random effects — alpha/gamma non-identifiable

</div>

<div class="p-2 rounded-lg bg-red-50 text-xs border border-red-300">

### <span class="text-red-600">v4–v5</span>
**Wrong turns**

v4: NegBin, dropped GP — supervisor rejected (cannot model spatial residuals)

v5: 3D product kernel GP — 100% bad HMC (250 latent vars)

</div>

<div class="p-2 rounded-lg bg-green-50 text-xs border border-green-400">

### <span class="text-green-600">v6 — Current</span>
**Spatial GP + AR(1)**

Matches epiwave.mapping exactly

I\* = rate (m*a*b*z), population in Poisson

Dual likelihood, <1% bad transitions

**Baseline run completed** — convergence issues diagnosed

</div>

</div>

<div class="mt-3 p-3 bg-blue-50 text-sm rounded-lg border border-blue-200">

**Key lessons:** Spatial GP + AR(1) from epiwave.mapping (fixes v1–v5 latent space), dual likelihood (addresses alpha/gamma identifiability), I\* as rate not count, WITH offset vs I\*=0 comparison as the core validation. Current work: diagnosing why alpha/gamma ridge persists and phi is underestimated.

</div>

<!--
Six iterations to get here. Versions 1-3 tried various GP approaches that failed computationally. Version 4 dropped the GP for a Negative Binomial approach — supervisor rejected it because it cannot model spatial residuals; it assumes I-star perfectly explains spatial variation. Version 5 used a 3D product kernel GP which was mathematically correct but computationally intractable — 100% bad HMC transitions because the latent space was too large.

Version 6 — our current design — follows epiwave.mapping directly. Spatial GP innovations with AR(1) temporal correlation keeps the latent space at n_sites. We've completed a baseline simulation-estimation run. The pipeline works — less than 1% bad HMC transitions — but we have diagnosed specific convergence issues with alpha-gamma identifiability and the spatial lengthscale phi that we'll discuss next.
-->
