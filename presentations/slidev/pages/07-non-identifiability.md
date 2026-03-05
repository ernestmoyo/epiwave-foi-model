---
layout: default
---

# Critical Discovery: Parameter Non-Identifiability

$\alpha$ and $\gamma$ trade off freely — only their **product** is estimable from count data

<div class="grid grid-cols-2 gap-4 mt-3">

<div class="p-3 rounded-lg bg-red-50 border border-red-300 text-sm">

### The Problem

Original model: $\mu = \gamma \cdot e^\alpha \cdot I^*$

For any constant $k$:

$$(\alpha + k,\; \gamma \cdot e^{-k}) \;\rightarrow\; \text{identical likelihood}$$

<div class="mt-3 text-sm">

| Chain | $\alpha$ | $\gamma$ | Product |
|:-----:|:--------:|:--------:|:-------:|
| 1 | -0.44 | 0.156 | **0.1003** |
| 2 | -0.50 | 0.165 | **0.1003** |

Chains wander the **ridge** of constant product

</div>
</div>

<div class="p-3 rounded-lg bg-green-50 border border-green-400 text-sm">

### The Fix

Merge into single identified parameter:

$$\lambda = \log(\gamma) + \alpha$$

$$e^\lambda = \text{reporting rate}$$

<div class="mt-2 text-sm">

**Result:**
- $\hat{R} < 1.05$
- ESS $> 1000$
- Clean convergence

</div>

<div class="mt-2 p-2 bg-white rounded border text-sm">

Also fixed: NegBin parameterisation **inverted** — `prob = mu/(mu+size)` should be `size/(size+mu)`. One character, hours of debugging.

</div>

</div>

</div>

<!--
This was a critical discovery during development. The original model had both an intercept alpha and a reporting rate gamma. These are fundamentally non-identifiable — for any constant k, you can shift alpha up by k and multiply gamma by e-to-the-minus-k and get exactly the same likelihood.

In practice, the MCMC chains would wander along this ridge. Chain 1 would find alpha equals negative 0.44 with gamma 0.156, Chain 2 would find alpha negative 0.5 with gamma 0.165 — but the product was always 0.1003, which is the true reporting rate of 0.1.

The fix was simple once we understood the problem — merge them into a single parameter lambda equals log(gamma) plus alpha. Now there's one identified parameter and convergence is clean.

I also want to mention a bug that cost hours to find — the Negative Binomial parameterisation was inverted. In greta, prob = size/(size+mu), not mu/(mu+size). One character difference. Every proposed likelihood value was near impossible, giving 0% acceptance. This is the kind of detail that matters enormously in practice.
-->
