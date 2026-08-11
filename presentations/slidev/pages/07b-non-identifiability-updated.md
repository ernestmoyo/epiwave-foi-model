---
layout: default
---

# Non-Identifiability: Resolved <span class="text-sm align-middle px-2 py-0.5 rounded-full bg-emerald-100 text-emerald-700 font-bold">UPDATED · Aug 2026</span>

<div class="grid grid-cols-3 gap-3 mt-2">

<div class="p-2 rounded-lg bg-red-50 border border-red-300 text-xs">

### The Problem (case-only)

cases ~ Poisson(γ · exp(α) · I\* · N)

(α + k, γ · exp(−k)) → **identical likelihood**

Individual α, γ not identifiable.

</div>

<div class="p-2 rounded-lg bg-green-50 border border-green-300 text-xs">

### The Fix — now working

Two pieces, both in since March:

1. **Dual likelihood** — prevalence (Binomial) informs α independently of γ.
2. **Prevalence depends on I** (GP-adjusted), not mechanistic X — so it actually constrains α.

**50-rep result:** α coverage **80%** WITH offset, **0%** for I\*=0.

</div>

<div class="p-2 rounded-lg bg-blue-50 border border-blue-300 text-xs">

### March question, answered

*"Is 30% survey coverage enough?"*

Yes — **once prevalence depends on I.** The offset is what anchors α's *level*; without it (I\*=0) α is unrecoverable even with the same surveys.

</div>

</div>

<div class="mt-3 p-2 bg-amber-50 border-l-4 border-amber-500 text-xs rounded-lg">

**What changed &amp; why:** In March this slide read *"ridge persists — dual likelihood has not resolved it."* The missing piece was that the prevalence likelihood was using mechanistic X, which is **independent of α** — so it could not break the ridge. Nick's fix (prevalence → I, the GP-adjusted incidence) makes the prevalence data depend on α, and the 50-rep coverage confirms the ridge is now broken.

</div>

<!--
This is the slide that changed the most in substance. In March I told you the alpha–gamma ridge persisted despite the dual likelihood, and I asked whether we simply needed more prevalence surveys. That turned out to be the wrong diagnosis. The real problem, which Nick identified, was that the prevalence likelihood was written against the mechanistic prevalence X — and X does not depend on alpha, so no amount of prevalence data could break the ridge. The fix is to make prevalence depend on I, the GP-adjusted latent incidence, which does carry alpha. With that change, the 50-replicate study shows alpha recovered at 80% coverage with the offset and zero percent without. So the answer to my March question — is 30% survey coverage enough — is yes, but only once the prevalence likelihood is wired to I. The offset then anchors the level of alpha.
-->
