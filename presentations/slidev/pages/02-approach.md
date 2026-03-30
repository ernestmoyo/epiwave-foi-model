---
layout: default
---

# Our Approach: Two-Stage Framework

Pre-compute transmission dynamics, then calibrate to data with a GP for spatial residuals

<div class="grid grid-cols-2 gap-6 mt-4">

<div class="p-3 rounded-lg bg-blue-50 border-l-4 border-blue-500 text-sm">

## Stage 1: Mechanistic Engine
**Ross-Macdonald ODE System**

Fixed entomological parameters from Vector Atlas:
- $m$ — mosquito density ratio
- $a$ — biting rate
- $g$ — mortality rate

<div class="mt-2 font-mono text-xs bg-white p-2 rounded">
Output: I*(s,t) = m·a·b·z — incidence rate
</div>

</div>

<div class="p-3 rounded-lg bg-green-50 border-l-4 border-green-500 text-sm">

## Stage 2: GP + Dual Likelihood
**Gaussian Process residuals + Poisson/Binomial**

Uses $I^*(s,t)$ as structural log-offset:

$$\log(I_{s,t}) = \alpha + \log(I^*_{s,t}) + \varepsilon_{s,t}$$

$$\varepsilon: \text{Spatial GP (Matérn 5/2) + AR(1) temporal}$$

<div class="mt-2 font-mono text-xs bg-white p-2 rounded">
5 parameters to infer via MCMC
</div>

</div>

</div>

<div class="mt-3 p-2 bg-gray-50 rounded text-center text-xs">

**Data flow:** Vector Atlas (m, a, g) → <code>apply_interventions()</code> → ODE solver → I*(s,t) → GP + dual likelihood → posterior predictions (5 params via HMC)

</div>

<!--
Our approach splits the problem into two stages. Stage 1 is the mechanistic engine — we take fixed entomological parameters from Vector Atlas or temperature-dependent models, and solve the Ross-Macdonald ODE system ONCE per spatial unit. This gives us I-star, the mechanistic prediction of expected incidence based purely on transmission dynamics.

Stage 2 then takes this mechanistic prediction as a structural offset in a GP model. The Gaussian Process captures spatially-correlated departures between the mechanistic prediction and what we observe in the data. We use a dual likelihood — Poisson for case counts and Binomial for prevalence surveys — which makes alpha (intercept) and gamma (reporting rate) individually identifiable. This follows the epiwave.mapping framework.
-->
