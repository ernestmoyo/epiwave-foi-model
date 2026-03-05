---
layout: two-cols
---

# Stage 1: Ross-Macdonald ODE System

Vector-informed transmission dynamics with intervention effects

<div class="mt-2 text-sm">

**Coupled ODE System:**

$$\frac{dx}{dt} = m \cdot a \cdot b \cdot z \cdot (1-x) - r \cdot x$$

$$\frac{dz}{dt} = a \cdot c \cdot x \cdot (1-z) - g \cdot z$$

**Force of Infection:**

$$\text{FOI}(t,s) = m(t,s) \cdot a(t,s) \cdot b \cdot z(t,s)$$

**Mechanistic Prediction:**

$$I^*(t,s) = \text{FOI}(t,s) \times N(s)$$

</div>

::right::

<div class="ml-4 mt-4 text-xs">

| Parameter | Description | Source |
|:---------:|:------------|:-------|
| $x$ | Human infection prevalence | ODE state |
| $z$ | Mosquito infection prevalence | ODE state |
| $m$ | Mosquito-to-human ratio | **Vector Atlas** |
| $a$ | Biting rate (bites/day) | **Vector Atlas** |
| $b$ | Prob: mosquito→human | Literature |
| $c$ | Prob: human→mosquito | Literature |
| $g$ | Mosquito mortality rate | **Vector Atlas** |
| $r$ | Human recovery rate | Literature |

<div class="mt-2 p-2 bg-amber-50 rounded-lg text-xs border border-amber-200">

All vector parameters **fixed** from Vector Atlas / entomological surveys — **no MCMC required**

</div>

</div>

<!--
Here's the mathematical detail of Stage 1. We have the classic Ross-Macdonald coupled ODE system. The dx/dt equation describes human infection dynamics — new infections from mosquito bites minus recovery. The dz/dt equation describes mosquito infection dynamics — new infections from biting infected humans minus mosquito death.

The key insight is that all the entomological parameters — m, a, and g — are FIXED from external data sources like Vector Atlas. We don't need to infer them. We solve these ODEs once per spatial unit to get z(t,s), then compute the mechanistic prediction I-star as the Force of Infection times population size.

The parameters b, c, and r are set from published literature values. Only the vector parameters that vary spatially come from Vector Atlas estimates.
-->
