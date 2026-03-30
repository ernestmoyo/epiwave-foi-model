---
layout: two-cols
---

# Stage 1: Ross-Macdonald ODE

Vector-informed transmission dynamics with intervention effects

<div class="mt-1 text-xs">

**Coupled ODE system:**

$$\frac{dx}{dt} = m \cdot a \cdot b \cdot z \cdot (1-x) - r \cdot x$$

$$\frac{dz}{dt} = a \cdot c \cdot x \cdot (1-z) - g \cdot z$$

**Mechanistic prediction (infection incidence RATE):**

$$I^*(t,s) = m(t,s) \cdot a(t,s) \cdot b \cdot z(t,s)$$

<div class="mt-2 p-2 bg-amber-50 rounded border border-amber-200 text-xs">
All vector parameters <strong>fixed</strong> from Vector Atlas / entomological surveys — <strong>no MCMC required</strong>. ODE solved <strong>once</strong> per site.
</div>

</div>

::right::

<div class="ml-2 mt-6">
<img src="/images/anim_ode_dynamics.gif" class="h-32 rounded shadow-lg mb-2" />
</div>

<div class="ml-2 text-xs">

| Parameter | Description | Source |
|:---------:|:------------|:-------|
| $x$ | Human prevalence | ODE state |
| $z$ | Mosquito prevalence | ODE state |
| $m$ | Mosquito-to-human ratio | **Vector Atlas** |
| $a$ | Biting rate (bites/day) | **Vector Atlas** |
| $b$ | Prob: mosquito to human | Literature (0.8) |
| $c$ | Prob: human to mosquito | Literature (0.8) |
| $g$ | Mosquito mortality rate | **Vector Atlas** |
| $r$ | Human recovery rate | Literature (1/7) |

<div class="mt-1 p-1 bg-blue-50 rounded text-xs border border-blue-200">
<strong>Interventions</strong> (ITN/IRS) adjust m, a, g <strong>upstream</strong> of ODE via <code>apply_interventions()</code>
</div>

</div>

<!--
Here's the mathematical detail of Stage 1. We have the classic Ross-Macdonald coupled ODE system. The dx/dt equation describes human infection dynamics — new infections from mosquito bites minus recovery. The dz/dt equation describes mosquito infection dynamics — new infections from biting infected humans minus mosquito death.

The key insight is that all the entomological parameters — m, a, and g — are FIXED from external data sources like Vector Atlas. We don't need to infer them. We solve these ODEs once per spatial unit to get z(t,s), then compute the mechanistic prediction I-star as the infection incidence rate: m times a times b times z. I-star is a rate, not a count — population enters the Poisson likelihood separately.

Interventions like ITN and IRS adjust the entomological parameters upstream, before the ODE is solved. So the ODE automatically reflects intervention effects.
-->
