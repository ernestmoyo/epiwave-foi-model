---
layout: default
---

# Stage 2 Evolution: GP to Negative Binomial

Three iterations of the residual structure, each solving the previous failure mode

<div class="grid grid-cols-4 gap-2 mt-4">

<div class="p-2 rounded-lg bg-red-50 text-xs border border-red-300">

### <span class="text-red-600">FAILED</span> v1
**Joint GP**

490 $\times$ 490 covariance matrix

Cholesky decomposition singular

<div class="mt-2 font-bold text-red-600 text-center">
0% HMC acceptance
</div>

</div>

<div class="p-2 rounded-lg bg-red-50 text-xs border border-red-300">

### <span class="text-red-600">FAILED</span> v2
**Separable GP**

10 $\times$ 10 + 49 $\times$ 49

Fixed matrix size but 59-dim latent space

<div class="mt-2 font-bold text-red-600 text-center">
HMC couldn't tune
</div>

</div>

<div class="p-2 rounded-lg bg-red-50 text-xs border border-red-300">

### <span class="text-red-600">FAILED</span> v3
**Hierarchical RE**

Normal random effects, non-centred param.

Poisson + 59 params

<div class="mt-2 font-bold text-red-600 text-center">
Extreme gradients
</div>

</div>

<div class="p-2 rounded-lg bg-green-50 text-xs border border-green-400">

### <span class="text-green-600">SUCCESS</span> v4
**NegBin**

2 parameters only

Overdispersion absorbs residual variation

<div class="mt-2 font-bold text-green-600 text-center">
90%+ acceptance
</div>

</div>

</div>

<div class="mt-3 p-3 bg-blue-50 text-sm rounded-lg border border-blue-200">

**Key Insight:** The mechanistic offset already encodes spatial-temporal structure — the GP was re-learning what the ODE already provides

</div>

<!--
This slide tells the story of how we arrived at the Negative Binomial formulation. It wasn't the first thing we tried — we went through three failed iterations.

Version 1 tried a joint space-time GP with a full 490 by 490 covariance matrix for 10 sites times 49 months. The Cholesky decomposition was near-singular, giving 0% HMC acceptance regardless of tuning.

Version 2 used a separable GP — a Kronecker product of a 10x10 spatial kernel and a 49x49 temporal kernel. This fixed the matrix size issue but created a 59-dimensional latent space that HMC simply couldn't navigate efficiently.

Version 3 tried hierarchical random effects with non-centred parameterisation. Still 59 parameters, and the Poisson likelihood with extreme count data created gradient issues.

The breakthrough was Version 4 — realising that the mechanistic offset already captures the spatial-temporal structure. The GP was trying to re-learn what the ODE already provides. A simple Negative Binomial with just 2 parameters gives 90%+ HMC acceptance.
-->
