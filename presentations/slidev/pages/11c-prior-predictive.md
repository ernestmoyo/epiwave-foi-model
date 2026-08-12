---
layout: default
---

# Prior Predictive Check: A Right Diagnosis, and a Deeper Problem <span class="text-sm align-middle px-2 py-0.5 rounded-full bg-emerald-100 text-emerald-700 font-bold">NEW · Aug 2026</span>

I ran Punam's prior predictive check. It correctly identified that the lengthscale φ was unidentifiable — but fixing φ made the mixing dramatically **worse**, which is how I found the real cause.

<div class="grid grid-cols-5 gap-3 mt-2">

<div class="col-span-3">
<img src="/images/ppc_identifiable_range.png" class="rounded shadow" />
</div>

<div class="col-span-2 text-xs flex flex-col gap-2">

<div class="p-2 rounded-lg bg-gray-50 border border-gray-300">
<b>The test.</b> For each draw of φ, compute the spatial correlation between the two <i>farthest</i> sites. Near 1 means every site is correlated with every other, the likelihood is flat in φ, and no sampling budget can identify it.
</div>

<div class="p-2 rounded-lg bg-red-50 border border-red-300">
<b>What it found.</b> φ = 3 implies <b>90%</b> correlation between the farthest sites — the truth sat inside the unidentifiable region, and 36% of the prior sat there too. Identifiable structure needs φ below about 1.0.
</div>

<div class="p-2 rounded-lg bg-amber-50 border border-amber-400">
<b>What happened when I fixed it.</b> Setting φ = 0.5 made φ identifiable and broke sampling outright — <b>R-hat 30–80, ESS ≈ 7</b> on the Workbench, across both dense and sparse GP configs. Far worse than the problem I was fixing.
</div>

<div class="p-2 rounded-lg bg-blue-50 border border-blue-300">
<b>The real cause.</b> σ² and θ trade off against each other: the AR(1) marginal variance is roughly <b>σ²/(1−θ²)</b>, so many pairs give the same field variance — and θ has a flat prior pinning nothing. A variance ridge, not a φ problem.
</div>

</div>

</div>

<div class="mt-2 p-2 bg-amber-50 border-l-4 border-amber-500 text-xs rounded-lg">

<b>Where this leaves things.</b> I have reverted to the working configuration (φ = 3) so the 50-replicate study still reproduces the valid result — α recovered in 76% of runs with the offset, 0% without. That result never depended on φ; it turns on the mechanistic offset and the dual likelihood. The GP hyperparameters remain unrecovered, and the fix is now specific: reparameterise the σ²/θ ridge and put a real prior on θ, rather than chasing φ. The check also flagged something I have not fixed — the priors imply a median prevalence of 0.97, so the incidence scale is still too high.

</div>

<!--
This slide is new, and it is the one where I want to be careful about what I am claiming. Punam suggested the prior predictive check as the way into the GP mixing problem, and I ran it. The test is on the right — for each draw of the lengthscale phi, compute the spatial correlation between the two farthest sites in the grid. If that number is near one, every site is effectively correlated with every other, the likelihood is flat in phi, and no sampling budget can recover it. What it found is that my simulation used phi equals three, which implies ninety percent correlation between the farthest sites. So I had placed the truth itself inside the unidentifiable region, and thirty-six percent of the prior was in there as well. That diagnosis was right. But here is the part worth reporting honestly — when I fixed it, setting phi to zero point five so that it became identifiable, the sampling did not improve, it collapsed. R-hat of thirty to eighty, effective sample size around seven, on the Workbench, across both the dense and the sparse GP configurations. Very much worse than the problem I was trying to fix. Chasing that down on the Workbench is what found the actual cause, and it is not phi at all. Sigma-squared and theta trade off against each other, because the marginal variance of the AR-one process is roughly sigma-squared over one minus theta-squared. Many combinations of the two give the same field variance, and theta has a flat prior, so nothing pins it down. That is a variance ridge, and it is a different problem from the lengthscale. So I have reverted to the working configuration with phi equal to three, which keeps the model sampleable and means the fifty-replicate result still stands — alpha recovered seventy-six percent of the time with the offset, never without it. That result never depended on phi. What I have gained is a specific target: reparameterise the sigma-squared theta ridge and put a real prior on theta, instead of chasing the lengthscale. And the check flagged one thing I have not fixed, which is that the priors imply a median prevalence of ninety-seven percent, so the incidence scale is still too high.
-->
