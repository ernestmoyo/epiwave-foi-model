---
layout: default
---

# MCMC Diagnostics: Narrowed to One Problem <span class="text-sm align-middle px-2 py-0.5 rounded-full bg-emerald-100 text-emerald-700 font-bold">UPDATED · Aug 2026</span>

<div class="grid grid-cols-5 gap-3 mt-2">

<div class="col-span-3">
<img src="/images/rhat_50rep.png" class="rounded shadow" />
</div>

<div class="col-span-2 text-xs flex flex-col gap-2">

<div class="p-2 rounded-lg bg-gray-50 border border-gray-300">
<b>What R-hat means.</b> R-hat compares the separate MCMC chains. When the chains agree it sits near 1; a value above about 1.1 means the chains have not settled on the same answer, so those estimates are not yet trustworthy.
</div>

<div class="p-2 rounded-lg bg-green-50 border border-green-300">
<b>What now mixes.</b> After running 4 chains of 2000 iterations (Punam), α and γ have R-hat between 1.02 and 1.10, with 500–790 effective samples. Their mixing is fine.
</div>

<div class="p-2 rounded-lg bg-red-50 border border-red-300">
<b>What does not.</b> The GP hyperparameters φ, σ² and θ still have R-hat between 2.2 and 3.7 across all 50 datasets. Because this holds over 50 datasets, it is a real mixing problem, not a one-run fluke.
</div>

</div>

</div>

<div class="mt-2 p-2 bg-amber-50 border-l-4 border-amber-500 text-xs rounded-lg">

<b>Why, and the plan.</b> The cause is in the simulation: with the spatial lengthscale φ set to 3 on coordinates rescaled to a 0–1 square, every pair of sites is more than 86% correlated, so the data cannot tell one lengthscale from another. Increasing iterations fixed α and γ but cannot fix this. The next step is the prior predictive check (Punam): draw φ, σ² and θ from their priors, simulate data, and check the priors put weight on sensible spatial patterns before fitting. I also plan a penalised-complexity prior on φ (Simpson et al.) to keep it from wandering.

</div>

<!--
In March I put three questions on this slide about convergence. Here is where they stand. R-hat is a convergence check — it compares the separate MCMC chains, sits near one when they agree, and above about one point one when they have not settled on the same answer. After running four chains of two thousand iterations, which was Punam's suggestion, alpha and gamma now mix well — R-hat near one, healthy effective sample sizes. But the GP hyperparameters, phi, sigma-squared and theta, still do not, with R-hat between two point two and three point seven, and they do so across all fifty datasets, which tells me it is a real mixing problem rather than a one-run fluke. The cause is in the simulation design: with phi set to three on coordinates rescaled to a zero-to-one square, every pair of sites is more than eighty-six percent correlated, so the data genuinely cannot distinguish one lengthscale from another. More iterations cannot fix that. The plan is Punam's prior predictive check, and a penalised-complexity prior on phi.
-->
