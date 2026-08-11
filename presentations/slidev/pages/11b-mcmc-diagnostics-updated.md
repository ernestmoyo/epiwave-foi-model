---
layout: default
---

# MCMC Diagnostics: Narrowed to One Problem <span class="text-sm align-middle px-2 py-0.5 rounded-full bg-emerald-100 text-emerald-700 font-bold">UPDATED · Aug 2026</span>

<div class="grid grid-cols-5 gap-3 mt-2">

<div class="col-span-3">
<img src="/images/rhat_50rep.png" class="rounded shadow" />
</div>

<div class="col-span-2 text-xs flex flex-col gap-2">

<div class="p-2 rounded-lg bg-green-50 border border-green-300">
<div class="font-bold text-green-700">March questions — answered</div>
Ran <b>2000 × 4 chains</b> (Punam's advice). <b>α, γ now mix</b>: R̂ 1.02–1.10, ESS 500–790.
</div>

<div class="p-2 rounded-lg bg-red-50 border border-red-300">
<div class="font-bold text-red-700">The one remaining issue</div>
GP hyperparams φ, σ², θ still stuck (R̂ 2.2–3.7, ESS 18–34). 50 reps → <b>structural, not sample-size</b>.
</div>

<div class="p-2 rounded-lg bg-gray-50 border border-gray-200">
<div class="font-bold">Root cause (stands from March)</div>
φ = 3 on [0,1] coords → &gt;86% correlation between all sites; data can't resolve the lengthscale.
</div>

</div>

</div>

<div class="mt-2 p-2 bg-amber-50 border-l-4 border-amber-500 text-xs rounded-lg">

**What changed &amp; why — and the plan:** More iterations resolved α/γ mixing but confirmed the GP-hyperparameter problem is **structural** (50 reps rules out noise). Next step is exactly **Punam's other suggestion — prior predictive checks**: draw φ, σ², θ from priors, simulate data, check the priors put mass on plausible spatial structure before fitting. Likely action: a **PC prior on φ** (Simpson et al.) to stop it wandering.

</div>

<!--
In March I asked whether to increase samples, adjust the true phi, and add more prevalence data. Here are the answers, and the plot makes them visual. This is median R-hat across fifty replicates — the dashed red line at 1.1 is the convergence threshold. Alpha and gamma, on the left, sit right at the line: increasing to four chains at two thousand iterations fixed their mixing. But the GP hyperparameters — sigma-squared, phi, theta — tower above at R-hat two to three point seven, and they do so in both models. Fifty replicates matters: it tells us this is a structural mixing problem, not a one-run fluke. The root cause is the one I diagnosed in March — phi equal to three on normalised coordinates makes every site more than 86 percent correlated, so the data genuinely cannot pin down the lengthscale. The plan is Punam's other suggestion, prior predictive checks, plus a penalised-complexity prior on phi. That directly closes her action item.
-->
