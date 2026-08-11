---
layout: default
---

# Next Steps: Updated <span class="text-sm align-middle px-2 py-0.5 rounded-full bg-emerald-100 text-emerald-700 font-bold">UPDATED · Aug 2026</span>

Status now: **offset validated across 50 reps** (α–γ ridge broken); remaining work isolated to GP-hyperparameter mixing + real data

<div class="grid grid-cols-2 gap-3 mt-3">

<div class="p-3 rounded-lg bg-blue-50 border-l-4 border-blue-500 text-sm">

### Priority plan

1. **Prior predictive checks** *(Punam's suggestion)* — the #1 item; targets the GP-hyperparameter mixing directly. PC prior on φ.
2. **Expand `apply_interventions()`** — biting-rate × human behaviour (Nick: the PhD contribution).
3. **Real data** — Mozambique / Angola pipeline; apply the validated framework.
4. **Framework paper (Obj. 2)** — now with a positive 50-rep result to anchor it.

</div>

<div class="p-3 rounded-lg bg-green-50 border-l-4 border-green-500 text-sm">

### March questions → resolved

- ✓ **α–γ ridge** — broken (prevalence-on-I fix).
- ✓ **Sampling budget** — 2000 × 4 chains; α/γ mix.
- ◐ **Simulation φ** — confirmed unidentifiable on [0,1] coords; prior-side fix planned.
- → **Real-data readiness** — Mozambique data on the MAP Workbench; scoping next.

</div>

</div>

<div class="mt-3 p-2 bg-amber-50 border-l-4 border-amber-500 text-xs rounded-lg">

**What changed &amp; why:** In March the plan was "address convergence, unclear how." It's now specific: α/γ are solved, so the entire convergence effort is focused on **one thing** — GP-hyperparameter mixing — with a named method (**prior predictive checks**, Punam) and a named prior (PC prior on φ). The framework itself is validated, so real-data application and the paper move up the queue.

</div>

<!--
The plan has sharpened since March. Back then the top line was "address convergence" without a clear method. Now alpha and gamma are solved, so the whole convergence effort narrows to a single target — the GP hyperparameters — and I have a specific method for it: Punam's prior predictive checks, plus a penalised-complexity prior on the lengthscale. Because the framework is now validated in simulation, the next real milestones are applying it to Mozambique or Angola data on the MAP Workbench, and writing the Objective 2 framework paper, which now has a genuinely positive fifty-replicate result to build on. The intervention work — modelling how biting rate responds to human behaviour — is the piece Nick flagged as the core PhD contribution and comes next after the mixing fix.
-->
