---
layout: default
---

# Our Approach: Two-Stage Framework

Pre-compute transmission dynamics, then calibrate to data — achieving **17x speedup** over joint inference

<div class="grid grid-cols-2 gap-6 mt-4">

<div class="p-3 rounded-lg bg-blue-50 border-l-4 border-blue-500 text-sm">

## Stage 1: Mechanistic Engine
**Ross-Macdonald ODE System**

Fixed entomological parameters from Vector Atlas:
- $m$ — mosquito density ratio
- $a$ — biting rate
- $g$ — mortality rate

<div class="mt-2 font-mono text-xs bg-white p-2 rounded">
Output: I*(s,t) — mechanistic prediction
</div>

</div>

<div class="p-3 rounded-lg bg-green-50 border-l-4 border-green-500 text-sm">

## Stage 2: Bayesian Calibration
**Negative Binomial Likelihood**

Uses $I^*(s,t)$ as structural log-offset:

$$\log(\mu) = \lambda + \log(I^*)$$

$$C \sim \text{NegBin}(\text{size}, p)$$

<div class="mt-2 font-mono text-xs bg-white p-2 rounded">
Only 2 parameters to infer via MCMC
</div>

</div>

</div>

<div class="mt-4 text-center p-2 bg-gray-100 rounded-lg">

**Key Advantage:** Mechanistic model is pre-computed and fixed — no MCMC over ODE parameters

</div>

<!--
Our approach splits the problem into two stages. Stage 1 is the mechanistic engine — we take fixed entomological parameters from Vector Atlas or temperature-dependent models, and solve the Ross-Macdonald ODE system ONCE per spatial unit. This gives us I-star, the mechanistic prediction of expected incidence based purely on transmission dynamics.

Stage 2 then takes this mechanistic prediction as a structural offset in a Bayesian calibration model. We use a Negative Binomial likelihood where the log-mean includes log(I-star) as a fixed offset. This means MCMC only needs to infer 2 parameters — the log reporting rate and the overdispersion — instead of the 20+ parameters required for joint inference. This is where the 17x speedup comes from.
-->
