---
layout: default
---

# Prior Predictive Check: The Mixing Problem Diagnosed <span class="text-sm align-middle px-2 py-0.5 rounded-full bg-emerald-100 text-emerald-700 font-bold">NEW · Aug 2026</span>

I ran Punam's prior predictive check — draw the parameters from their priors, push them through the model, and look at the data they imply, before fitting anything. It found the cause of the GP mixing problem, and it was in my simulation, not in the model.

<div class="grid grid-cols-5 gap-3 mt-2">

<div class="col-span-3">
<img src="/images/ppc_before.png" class="rounded shadow" />
</div>

<div class="col-span-2 text-xs flex flex-col gap-2">

<div class="p-2 rounded-lg bg-gray-50 border border-gray-300">
<b>The test.</b> For each draw of the lengthscale φ, compute the spatial correlation between the two <i>farthest</i> sites. If that number is near 1, every site is effectively correlated with every other, the likelihood is flat in φ, and no amount of sampling can identify it.
</div>

<div class="p-2 rounded-lg bg-red-50 border border-red-300">
<b>What it found.</b> My simulation used φ = 3, which implies <b>90%</b> correlation between the farthest sites — the truth was placed inside the unidentifiable region. The prior made it worse: <b>36%</b> of its mass sat there too.
</div>

<div class="p-2 rounded-lg bg-green-50 border border-green-300">
<b>The fix, now in the code.</b> The φ prior moved to a median of 0.5, putting <b>0%</b> of its mass in the unidentifiable region; the simulation truth moved from φ = 3 to φ = 0.5, which implies 11% correlation. On this grid, identifiable structure needs φ below about 1.0.
</div>

</div>

</div>

<div class="mt-2 p-2 bg-amber-50 border-l-4 border-amber-500 text-xs rounded-lg">

<b>What this changes, and what it does not.</b> It explains why more iterations fixed α and γ but did nothing for φ, σ² and θ — those three were unidentifiable by construction, so no sampling budget could have recovered them. It does <b>not</b> yet mean they are recovered: the 50-replicate study was run before this fix, so re-running it under φ = 0.5 is the next job. The α result is untouched by this, because it turns on the offset and the dual likelihood, not on the spatial kernel. The check also flagged one thing I have not fixed: the priors imply a median prevalence of 0.97, so the incidence scale is still too high.

</div>

<!--
This slide is new since the Melbourne discussion. Punam suggested the prior predictive check as the way into the GP mixing problem, and I have now run it. The idea is simple — before fitting anything, draw the parameters from their priors, push them through the model, and look at the data those priors imply. The specific test that mattered is on the right. For each draw of the lengthscale phi, I compute the spatial correlation between the two farthest sites in the grid. If that number is close to one, then every site is effectively correlated with every other site, the likelihood is flat in phi, and the data simply cannot tell one lengthscale from another. What it found is that my simulation used phi equals three, which implies ninety percent correlation between the farthest sites — so I had put the truth itself inside the unidentifiable region. And the prior compounded it, with thirty-six percent of its mass in the same place. Both are now fixed in the code: the prior median moved to zero point five, which puts none of its mass in the bad region, and the simulation truth moved to zero point five as well. The honest caveat is that this explains the mixing problem but does not yet solve it — the fifty-replicate study I showed you was run before this fix, so I need to re-run it before I can claim phi, sigma-squared and theta are recovered. Alpha is unaffected, because alpha depends on the offset and the dual likelihood rather than on the spatial kernel. And the check flagged one thing I have not yet fixed — the priors imply a median prevalence of ninety-seven percent, which is far too high, so the incidence scale needs tightening before this meets real data.
-->
