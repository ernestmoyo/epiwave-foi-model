---
layout: default
---

# Questions and Next Steps

Current status: pipeline works end-to-end, but convergence needs addressing before real data

<div class="grid grid-cols-2 gap-3 mt-3">

<div class="p-3 rounded-lg bg-amber-50 border-l-4 border-amber-500 text-sm">

### Questions for Supervisors

1. **Alpha-gamma ridge:** The dual likelihood has not fully broken the ridge at 30% survey coverage. How much prevalence data is needed? Is this a simulation design issue or a structural limitation?

2. **Simulation phi:** True phi=3.0 on normalised [0,1] coords means >86% spatial correlation between all sites. Should we lower phi or use raw coordinates so the GP has actual spatial variation to estimate?

3. **Sampling budget:** 1000 samples shows multimodal posteriors. Should we invest in 4000+ samples before diagnosing further, or address the priors first?

4. **Real data readiness:** What data sources should we prioritise for the study site — publicly available prevalence surveys, routine case data, or both?

</div>

<div class="p-3 rounded-lg bg-blue-50 border-l-4 border-blue-500 text-sm">

### Immediate Plan

1. **Address convergence** — increase samples, review prior-data consistency for phi and sigma2
2. **Sensitivity study** — vary survey_fraction (10%, 30%, 50%) to quantify how much prevalence data breaks the ridge
3. **Real data scoping** — identify available data for the study site

### Framework Paper (Obj. 2)

4. **Misspecification tests** — does the offset help even when I\* is wrong? (robustness check)
5. **Write-up** — methods, simulation results, WITH vs WITHOUT comparison


</div>

</div>

<!--
These are the questions I need guidance on. The key convergence issues are well-diagnosed — we know what's happening and have specific hypotheses about why. But I need input on which direction to prioritise.

The alpha-gamma ridge is the most important question. The dual likelihood should break it mathematically. But in our simulation with 147 prevalence surveys, the ridge persists. Is this because we don't have enough surveys, or is there something else going on?

The phi question is about simulation design — not the model itself. On normalised coordinates, phi=3 produces near-total spatial correlation, so there's almost no information in the data to estimate the lengthscale. We could either lower the true phi or use unnormalised coordinates.

On real data — we need to scope what's available for the study site and confirm data access permissions.
-->
