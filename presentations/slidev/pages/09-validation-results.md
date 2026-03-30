---
layout: default
---

# Simulation-Estimation: Current Results

Core comparison: **WITH I\* offset** vs **I\*=0 (standard geostatistical)** &nbsp;|&nbsp; 10 sites x 49 time steps, 1000 samples

<div class="grid grid-cols-2 gap-4 mt-4">

<div class="p-3 rounded-lg bg-green-50 text-center border border-green-300">
<div class="text-xl font-bold text-green-700">WITH I* Offset</div>
<div class="text-xs mt-1">log(I) = α + log(I*) + ε</div>
<div class="mt-2 text-xs">

| Param | Estimated | True | Status |
|:------|:---------:|:----:|:------:|
| $\alpha$ | -0.66 | 0.0 | biased |
| $\gamma$ | 0.15 | 0.10 | biased |
| $\sigma^2$ | 0.79 | 0.36 | overest. |
| $\phi$ | 1.53 | 3.0 | underest. |
| $\theta$ (AR1) | **0.73** | **0.75** | good |

</div>
<div class="mt-1 font-bold text-green-600">&lt;1% bad HMC transitions</div>
</div>

<div class="p-3 rounded-lg bg-blue-50 text-center border border-blue-300">
<div class="text-xl font-bold text-blue-700">I*=0 (Standard Geostatistical)</div>
<div class="text-xs mt-1">log(I) = α + ε</div>
<div class="mt-2 text-xs">

| Param | Estimated | True | Status |
|:------|:---------:|:----:|:------:|
| $\alpha$ | -2.16 | — | much lower |
| $\gamma$ | 0.12 | — | similar |
| $\sigma^2$ | **3.69** | 0.36 | 5x larger |
| $\phi$ | 1.24 | 3.0 | underest. |
| $\theta$ (AR1) | 0.61 | 0.75 | lower |

</div>
<div class="mt-1 font-bold text-amber-600">11% bad transitions (warmup)</div>
</div>

</div>

<div class="mt-3 p-2 bg-gray-50 text-xs rounded-lg">

**Key finding:** WITHOUT offset needs 5x more GP variance ($\sigma^2$ = 3.69 vs 0.79) — confirming that I\* provides spatial structure the GP would otherwise have to learn alone. AR(1) temporal correlation $\theta$ recovers well in the WITH model. **Open questions:** alpha and phi are not recovering true values — see diagnostics slide.

</div>

<!--
These are the actual numbers from our most recent run — 10 sites, 49 time steps, 1000 MCMC samples with 500 warmup, 2 chains.

The WITH model recovers AR(1) temporal correlation well — theta is 0.73 versus the true 0.75. But alpha is biased at -0.66 instead of 0, and gamma is 0.15 instead of 0.10. These are trading off along the alpha-gamma ridge we saw in the previous slide. Phi (spatial lengthscale) is estimated at 1.53 versus the true 3.0 — we'll discuss why in the diagnostics.

The WITHOUT model — I-star equals zero, which is the standard geostatistical approach — needs five times more GP variance. Sigma-squared jumps from 0.79 to 3.69 because without the mechanistic offset, the GP has to absorb all the spatial variation that I-star would have explained. It also has 11% bad HMC transitions during warmup compared to less than 1% for the WITH model.

This confirms the hypothesis: including the mechanistic prediction provides valuable spatial structure. But there are convergence issues to address before this is a convincing result.
-->
