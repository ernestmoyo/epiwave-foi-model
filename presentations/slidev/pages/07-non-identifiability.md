---
layout: default
---

# Parameter Non-Identifiability: Current Status

<div class="grid grid-cols-3 gap-3 mt-2">

<div class="p-2 rounded-lg bg-red-50 border border-red-300 text-xs">

### The Problem (case-only)

cases ~ Poisson(γ · exp(α) · I\* · N)

For any constant k:

(α + k, γ · exp(-k)) → **identical likelihood**

The product exp(α)·γ is identifiable, but individual parameters are not.

</div>

<div class="p-2 rounded-lg bg-amber-50 border border-amber-300 text-xs">

### The Fix: Dual Likelihood

Add **prevalence surveys** (Binomial):

$$Y_{s,t} \sim \text{Binomial}(N_{s,t},\; x_{s,t})$$

Prevalence is **unaffected by** γ, so it separately informs α.

**Current result (10 sites, 147 surveys):**
Ridge persists. **Question:** Is 30% survey coverage sufficient?

</div>

<div class="flex items-center justify-center">
<img src="/images/alpha_gamma_ridge.png" class="h-44 rounded shadow-lg" />
</div>

</div>

<div class="mt-2 text-xs text-center text-gray-600">
Dashed lines = true values (α=0, γ=0.1). The negative correlation ridge is visible — dual likelihood has not fully resolved it at this sample size.
</div>

<!--
This is a critical finding. With case data alone, alpha and gamma are fundamentally non-identifiable. The dual likelihood design — adding prevalence surveys via a Binomial likelihood — should break this ridge because prevalence data informs infection incidence independently of reporting rate.

However, our current simulation with 147 prevalence surveys across 490 site-time combinations still shows a clear ridge in the joint posterior. The question for supervisors is: is 30% survey coverage sufficient? With real data, how much prevalence data would we typically have, and is it enough to resolve this degeneracy? This is a legitimate research question — how much prevalence information does the dual likelihood need to work?
-->
