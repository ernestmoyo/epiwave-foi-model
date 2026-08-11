---
layout: default
---

# Simulation-Estimation: Updated Results <span class="text-sm align-middle px-2 py-0.5 rounded-full bg-emerald-100 text-emerald-700 font-bold">UPDATED · Aug 2026</span>

I simulated 50 datasets from the model, fitted both versions to each, and asked how often each fitted parameter recovered the true value used to simulate the data (50 replicates, 10 sites x 48 time steps, 2000 warmup / 2000 samples, 4 chains).

<div class="grid grid-cols-5 gap-3 mt-2">

<div class="col-span-3">
<img src="/images/coverage_50rep.png" class="rounded shadow" />
</div>

<div class="col-span-2 text-xs flex flex-col gap-2">

<div class="p-2 rounded-lg bg-gray-50 border border-gray-300">
<b>What the plot means.</b> α is the intercept on the log-incidence scale, so it sets the overall level of infection incidence. Coverage is the fraction of the 50 datasets in which the 95% credible interval for a parameter contained the true value. "80% coverage for α" means the interval contained the true α in 40 of the 50 runs.
</div>

<div class="p-2 rounded-lg bg-green-50 border border-green-300">
<b>Main result.</b> The model recovers α in 80% of runs when it uses the vector offset I*, and in none of the runs when the offset is removed (I* = 0). The vector information is what lets the model recover the overall level of incidence.
</div>

<div class="p-2 rounded-lg bg-blue-50 border border-blue-300">
<b>Reporting rate γ.</b> γ is recovered in about 80% of runs either way, because the prevalence surveys inform it directly regardless of the offset.
</div>

</div>

</div>

<div class="mt-2 p-2 bg-amber-50 border-l-4 border-amber-500 text-xs rounded-lg">

<b>What changed since March, and one caveat.</b> Since March I made three changes: I corrected the prevalence likelihood so it depends on the modelled incidence I rather than the mechanistic prevalence X (Nick); I centred α and ran 4 chains of 2000 iterations (Punam); and I scored 50 replicates so I report how often the truth is recovered, not one single run. Caveat raised in Melbourne: these runs assume 30% prevalence-survey coverage, which is very high for real settings (David) — the result needs re-testing with survey numbers taken from real prevalence data.

</div>

<!--
This replaces the single-run slide from March, and I have rewritten it so the terms are defined on the slide. Alpha is the intercept on the log-incidence scale — it sets the overall level of infection incidence. Coverage is the fraction of the fifty simulated datasets in which the ninety-five percent credible interval for a parameter contained the true value I used to simulate that dataset. So eighty percent coverage for alpha means that in forty of the fifty runs, the interval I estimated contained the true alpha. The plot shows this per parameter, green for the model with the vector offset and blue for the model without it. Alpha is recovered eighty percent of the time with the offset and zero percent without it — the vector information is what lets the model pin down the overall level of incidence. Gamma, the reporting rate, is recovered either way because the prevalence surveys inform it. The caveat David raised in Melbourne is that thirty percent prevalence-survey coverage is unrealistically high, so I need to re-run this with survey numbers taken from real prevalence data before claiming the result holds in practice.
-->
