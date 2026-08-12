---
layout: default
---

# Next Steps: Updated <span class="text-sm align-middle px-2 py-0.5 rounded-full bg-emerald-100 text-emerald-700 font-bold">UPDATED · Aug 2026</span>

The framework now recovers α and γ in simulation. The remaining work is the GP-hyperparameter mixing, a more realistic survey design, and real data.

<div class="grid grid-cols-2 gap-3 mt-3">

<div class="p-3 rounded-lg bg-blue-50 border-l-4 border-blue-500 text-sm">

### Plan, in order

1. Reparameterise the σ²/θ variance ridge and put a real prior on θ — the actual cause of the GP-hyperparameter mixing, found after the φ fix backfired.
2. Re-run the study with prevalence-survey numbers taken from real data, since 30% coverage is unrealistically high (David).
3. Tighten the α and σ² priors — the prior predictive check shows they imply a median prevalence of 0.97.
4. Expand `apply_interventions()` so the biting rate responds to human behaviour and intervention timing (Nick).
5. Apply the framework to real data (Mozambique or Angola) on the MAP Workbench.

</div>

<div class="p-3 rounded-lg bg-green-50 border-l-4 border-green-500 text-sm">

### March questions, where they stand

- The α–γ trade-off is resolved: prevalence now depends on the modelled incidence I (Nick).
- Sampling budget is settled: 4 chains of 2000 iterations, so α and γ mix.
- The lengthscale φ is confirmed unidentifiable at φ = 3 — but making it identifiable broke sampling, so the code is back on the working config and the real target is the σ²/θ ridge.
- The Objective 2 paper is updated and pushed, anchored on the 50-replicate result.
- Real-data readiness: the Mozambique data is on the MAP Workbench, and scoping is next.

</div>

</div>

<div class="mt-3 p-2 bg-amber-50 border-l-4 border-amber-500 text-xs rounded-lg">

<b>What changed since March.</b> In March the plan was simply "address convergence" without a method. Since then the convergence work has narrowed to one target, the GP hyperparameters, and Punam's prior predictive check has been run rather than merely named. It found that the lengthscale was unidentifiable — but correcting that broke sampling outright, which is how the real cause surfaced: a variance ridge between σ² and θ. The code is back on the working configuration, so the α result stands, and the remaining fix is now specific. The Melbourne discussion added a second item: the 30% survey coverage is too high to be realistic, so the identifiability result has to be re-tested against survey numbers drawn from real prevalence data (David).

</div>

<!--
The plan has sharpened since March. Back then the top line was address convergence, with no method. Now alpha and gamma are recovered, so the convergence work is focused on one target, the GP hyperparameters — and Punam's prior predictive check has been run rather than merely planned. It found that the lengthscale was unidentifiable, but when I corrected that the sampling collapsed, and chasing that down is what exposed the real cause: a variance ridge between sigma-squared and theta. I have reverted to the working configuration so the alpha result stands, and item one is now to reparameterise that ridge and put a prior on theta. The Melbourne discussion added a second item that I need to act on: thirty percent prevalence-survey coverage is unrealistically high, so before I present this as a practical result I have to re-run it with survey numbers taken from real prevalence data, which is David's point. Third is tightening the priors, since the same check showed they imply a median prevalence of ninety-seven percent. After that comes the intervention work Nick flagged as the core contribution, and then real data on the Workbench. The paper itself is updated and pushed, anchored on the fifty-replicate result.
-->
