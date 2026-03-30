---
layout: default
---

# MCMC Diagnostics: What's Working and What's Not

<div class="grid grid-cols-2 gap-3 mt-2">

<div>
<div class="text-xs font-bold mb-1">GP Hyperparameter Posteriors</div>
<img src="/images/gp_hyperparameter_posteriors.png" class="rounded shadow" />
<div class="text-xs mt-1 text-gray-600">Red dashed = true values. Multimodal posteriors for phi and sigma2 suggest insufficient mixing.</div>
</div>

<div>
<div class="text-xs font-bold mb-1">Posterior Predictive Check (Site 1)</div>
<img src="/images/posterior_predictive_check.png" class="rounded shadow" />
<div class="text-xs mt-1 text-gray-600">GP+offset captures the shape but underestimates peaks — consistent with alpha bias.</div>
</div>

</div>

<div class="mt-2 p-2 bg-amber-50 rounded-lg border border-amber-200 text-xs">

**Diagnosed issues (from walkthrough analysis):**

1. **phi = 3.0 on [0,1] normalised coords** produces >86% correlation between ALL site pairs — the data cannot distinguish phi=3 from phi=10. Prior (lognormal(0.5,0.5), median=1.65) pulls phi toward ~1.5.
2. **1000 samples / 500 warmup** may be insufficient — multimodal posteriors in phi and sigma2.
3. **147 prevalence surveys** across 490 site-time pairs may not break the alpha-gamma ridge fully.

**Questions for supervisors:** Should we increase samples (4000+)? Adjust the simulation's true phi to a value the data can actually distinguish? Increase survey coverage in the simulation?

</div>

<!--
These are the actual diagnostic plots from our run. The GP hyperparameter posteriors show multimodal distributions — phi concentrates around 1.2 to 1.5 but the true value is 3.0 (red dashed line). We diagnosed why: on normalised coordinates in a [0,1] square, phi equals 3 means minimum off-diagonal correlation of 86%. Every site is nearly perfectly correlated with every other site. The data simply cannot distinguish phi equals 3 from phi equals 5 or 10 — they all produce nearly identical correlation matrices.

The posterior predictive check shows the GP plus offset model captures the seasonal shape but consistently underestimates the peaks. This is consistent with the alpha bias — alpha is estimated at -0.66 instead of 0, which pulls all predictions down.

These are questions I need guidance on: Should I increase the MCMC samples? Should I adjust the simulation design so that the true phi is on a scale the data can actually resolve? And how much prevalence data do we need for the dual likelihood to break the alpha-gamma ridge?
-->
