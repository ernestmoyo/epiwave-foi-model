---
layout: default
---

# Simulation-Estimation: Updated Results <span class="text-sm align-middle px-2 py-0.5 rounded-full bg-emerald-100 text-emerald-700 font-bold">UPDATED · Aug 2026</span>

**50-replicate study** &nbsp;|&nbsp; 10 sites x 48 times &nbsp;|&nbsp; **2000 warmup / 2000 samples / 4 chains** &nbsp;|&nbsp; scored by **coverage**, not one point estimate

<div class="grid grid-cols-5 gap-3 mt-2">

<div class="col-span-3">
<img src="/images/coverage_50rep.png" class="rounded shadow" />
</div>

<div class="col-span-2 text-xs flex flex-col gap-2">

<div class="p-2 rounded-lg bg-green-50 border border-green-300">
<div class="font-bold text-green-700">α — the headline</div>
80% coverage <b>WITH</b> offset vs <b>0%</b> for I\*=0. The ridge is broken.
</div>

<div class="p-2 rounded-lg bg-blue-50 border border-blue-300">
<div class="font-bold text-blue-700">γ (reporting rate)</div>
80% / 78% — identified either way (prevalence pins it down).
</div>

<div class="p-2 rounded-lg bg-purple-50 border border-purple-200">
<div class="font-bold text-purple-700">Predictive edge</div>
Latent-surface RMSE <b>0.052</b> (WITH) vs <b>0.058</b> (I\*=0).
</div>

<div class="p-2 rounded-lg bg-red-50 border border-red-200">
<div class="font-bold text-red-700">Caveat → next slide</div>
GP hyperparams (σ², φ, θ) low coverage — not yet converged.
</div>

</div>

</div>

<div class="mt-2 p-2 bg-amber-50 border-l-4 border-amber-500 text-xs rounded-lg">

**What changed &amp; why since March:** (i) Nick's **prevalence-likelihood fix** (prevalence → GP-adjusted I, not mechanistic X); (ii) **centred α** + **4 chains × 2000** (Punam's iteration advice); (iii) **50 replicates** so we report *coverage*, not one possibly-lucky run. The single-run "α = −0.66, biased" became **α at 80% coverage WITH the offset vs 0% without.**

</div>

<!--
This replaces the single-run slide from March. The plot is the whole story: it shows 95% credible-interval coverage across fifty replicates for each parameter, green for the model with the mechanistic offset and blue for the standard geostatistical model with I-star set to zero. Look at alpha on the far left — 80% coverage with the offset, zero without it. That is the alpha–gamma ridge, broken. Gamma recovers either way because the prevalence likelihood pins it down. The offset also wins modestly on predictive RMSE. The honest caveat is the right-hand end of the plot: the GP hyperparameters sigma-squared, phi and theta have low coverage — they are not converging — and that is the subject of the next slide. Three things changed since March to get here: Nick's prevalence-on-I fix, centring alpha with four chains at two thousand iterations, and scoring fifty replicates so we report calibration rather than eyeballing one run.
-->
