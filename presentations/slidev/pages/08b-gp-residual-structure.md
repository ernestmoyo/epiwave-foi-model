---
layout: default
---

# True GP Residual Structure

The spatial GP + AR(1) creates the correlated epsilon that Stage 2 must recover

<div class="flex justify-center mt-2">
<img src="/images/gp_residual_heatmap.png" class="h-48 rounded shadow-lg" />
</div>

<div class="grid grid-cols-3 gap-3 mt-3 text-xs">

<div class="p-2 bg-blue-50 rounded border border-blue-200">

**Spatial correlation** (columns)
Nearby sites share similar colours at each time step — driven by Matern 5/2 kernel with $\phi = 3.0$.

</div>

<div class="p-2 bg-green-50 rounded border border-green-200">

**Temporal persistence** (rows)
Colours persist across adjacent months — driven by AR(1) with $\theta = 0.75$. Residuals are autocorrelated.

</div>

<div class="p-2 bg-amber-50 rounded border border-amber-200">

**What this represents**
epsilon captures where reality **differs** from I\*. Blue = true incidence higher than mechanistic prediction. Red = lower.

</div>

</div>

True parameters: σ = 0.6, ϕ = 3.0 (spatial), θ = 0.75 (AR1 temporal). 10 sites, 49 time steps.
{.text-xs .text-center .mt-2 .p-2 .bg-gray-50 .rounded}

<!--
This heatmap shows the true GP residuals that we simulated — this is what Stage 2 needs to recover. Each cell represents epsilon at a particular site and month. Blue means the true incidence is higher than the mechanistic prediction, red means lower.

You can see two patterns. First, spatial correlation — at any given time step, nearby sites tend to have similar colours. This comes from the Matern 5/2 kernel with lengthscale phi equals 3. Second, temporal persistence — the colours persist across adjacent months, which comes from the AR(1) temporal correlation with theta equals 0.75.

This structure is exactly what the framework aims to capture. The GP models how spatial variation in case incidence differs from the mechanistic prediction. Without this GP, the model would assume I-star perfectly explains spatial variation — which it clearly does not.
-->
