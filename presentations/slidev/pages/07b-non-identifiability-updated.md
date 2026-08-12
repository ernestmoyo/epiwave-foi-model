---
layout: default
---

# Non-Identifiability: Resolved <span class="text-sm align-middle px-2 py-0.5 rounded-full bg-emerald-100 text-emerald-700 font-bold">UPDATED · Aug 2026</span>

<div class="grid grid-cols-3 gap-3 mt-2">

<div class="p-2 rounded-lg bg-red-50 border border-red-300 text-xs">

### The problem (incidence from case counts only)

Case counts are modelled as a reporting fraction of infection incidence:

cases ~ Poisson(γ · exp(α) · I\* · N)

For any constant k, the pair (α + k, γ · exp(−k)) gives the same case counts. The product exp(α)·γ is identifiable from cases, but α and γ on their own are not.

</div>

<div class="p-2 rounded-lg bg-green-50 border border-green-300 text-xs">

### The fix, and why it now works

I add prevalence surveys as a second likelihood. Prevalence measures infection directly, so it does not depend on the reporting rate γ and can inform α on its own.

The correction (Nick): the prevalence likelihood must depend on the modelled incidence I, not the mechanistic prevalence X. X does not contain α, so on its own it could not break the tie.

</div>

<div class="p-2 rounded-lg bg-blue-50 border border-blue-300 text-xs">

### Result across 50 datasets

With the offset, the model recovers α in 76% of the simulated datasets. Without it (I\* = 0), α is recovered in none of them.

Caveat (Melbourne, David): this used 30% survey coverage, which is unrealistically high — to be re-tested with real survey numbers.

</div>

</div>

<div class="mt-3 p-2 bg-amber-50 border-l-4 border-amber-500 text-xs rounded-lg">

<b>What changed since March.</b> In March this slide said the ridge persisted despite the dual likelihood, and I asked whether more prevalence surveys would fix it. That was the wrong diagnosis. The prevalence likelihood was written against the mechanistic prevalence X, which does not contain α, so no amount of prevalence data could separate α from γ. Once prevalence is made to depend on the modelled incidence I (Nick), the surveys inform α directly and the 50-replicate study confirms α is now recovered.

</div>

<!--
This is the slide that changed the most in substance, and I have rewritten it to talk about infection incidence rather than only cases, to match the earlier approach slide. The model treats case counts as a reporting fraction of infection incidence. With case counts alone, alpha the incidence intercept and gamma the reporting rate trade off exactly — you can raise one and lower the other and get identical case counts — so only their product is identifiable. In March I told you the dual likelihood had not fixed this and asked whether we needed more surveys. That was the wrong diagnosis. The prevalence likelihood was written against X, the mechanistic prevalence, which does not contain alpha, so it could never separate alpha from gamma. Nick's correction is to make prevalence depend on I, the modelled incidence, which does contain alpha. With that change the fifty-replicate study recovers alpha seventy-six percent of the time with the offset and never without it. The caveat David raised is that thirty percent survey coverage is unrealistically high, so I need to re-test with realistic survey numbers.
-->
