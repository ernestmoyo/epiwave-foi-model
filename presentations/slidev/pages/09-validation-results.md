---
layout: default
---

# Simulation-Estimation Validation

10 sites $\times$ 49 months &nbsp;|&nbsp; 2000 MCMC samples &nbsp;|&nbsp; 2 chains &nbsp;|&nbsp; HMC via greta/TensorFlow

<div class="grid grid-cols-4 gap-3 mt-4">

<div class="p-3 rounded-lg bg-blue-50 text-center">
<div class="text-3xl font-bold text-blue-700">88%</div>
<div class="text-sm mt-1 font-semibold">RMSE Improvement</div>
<div class="text-xs mt-1 opacity-70">WITH offset vs WITHOUT</div>
</div>

<div class="p-3 rounded-lg bg-green-50 text-center">
<div class="text-3xl font-bold text-green-700">14.7</div>
<div class="text-sm mt-1 font-semibold">RMSE WITH offset</div>
<div class="text-xs mt-1 opacity-70">7.2% of mean observed</div>
</div>

<div class="p-3 rounded-lg bg-red-50 text-center">
<div class="text-3xl font-bold text-red-600">125.5</div>
<div class="text-sm mt-1 font-semibold">RMSE WITHOUT offset</div>
<div class="text-xs mt-1 opacity-70">61.3% of mean observed</div>
</div>

<div class="p-3 rounded-lg bg-purple-50 text-center">
<div class="text-3xl font-bold text-purple-700">0.1003</div>
<div class="text-sm mt-1 font-semibold">Recovered reporting rate</div>
<div class="text-xs mt-1 opacity-70">True = 0.1000</div>
</div>

</div>

<div class="mt-3 p-2 bg-gray-50 text-sm rounded-lg">

**Precision comparison:** WITH model 95% CI width = **0.0013** vs WITHOUT = **0.0133** — mechanistic offset provides **10x tighter** credible intervals

</div>

<div class="mt-2 text-xs">

| Metric | WITH offset | WITHOUT offset | Improvement |
|:-------|:-----------:|:--------------:|:-----------:|
| RMSE | 14.7 | 125.5 | **88%** |
| RMSE / mean | 7.2% | 61.3% | — |
| 95% CI width ($\lambda$) | 0.0013 | 0.0133 | **10x** |
| $\hat{R}$ (log_rate) | 1.004 | 1.002 | Both converged |
| Reporting rate recovery | 0.1003 | — | Near-zero bias |

</div>

<!--
Here are the key results from the simulation-estimation validation. The headline number is 88% RMSE improvement when using the mechanistic offset compared to the standard model without it.

The WITH offset model achieves an RMSE of 14.7, which is only 7.2% of the mean observed case count. Without the offset, the RMSE balloons to 125.5 — more than 61% of the mean.

Critically, the reporting rate is recovered almost exactly. The true value was 0.1, and we estimated 0.1003 — essentially zero bias. This validates the entire two-stage framework: Stage 1 correctly encodes transmission dynamics, and Stage 2 correctly recovers the scaling between mechanistic predictions and observed data.

The precision difference is also striking — the credible interval for the reporting rate is 10 times tighter with the mechanistic offset. This makes biological sense: when you tell the model the shape of the transmission curve, it only needs to estimate the scale, so uncertainty is dramatically reduced.
-->
