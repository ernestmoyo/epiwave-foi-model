---
layout: default
---

# The Challenge

Current malaria risk mapping treats vector biology and case data as **separate worlds**

<div class="grid grid-cols-3 gap-4 mt-4">

<div class="p-3 rounded-lg bg-blue-50 border border-blue-200 text-sm">

### Purely Statistical

Standard geostatistical models (MAP / Bhatt et al.) fit surfaces to case data **without mechanistic understanding** of transmission dynamics

</div>

<div class="p-3 rounded-lg bg-red-50 border border-red-200 text-sm">

### Purely Mechanistic

Full ODE-based approaches require MCMC over **all dynamic parameters** — computationally intractable at continental scale

</div>

<div class="p-3 rounded-lg bg-amber-50 border border-amber-200 text-sm">

### No Bridge

Vector surveillance data (biting rates, survival, resistance) from entomological studies remains **disconnected** from routine case reporting

</div>

</div>

<div class="mt-4 p-3 bg-gray-100 rounded-lg text-center italic">
"We need a framework that combines mechanistic understanding with statistical flexibility — without the computational cost"
</div>

<!--
The core challenge is threefold. First, purely statistical approaches like MAP don't incorporate the biological mechanisms of malaria transmission. Second, purely mechanistic approaches that try to infer all ODE parameters via MCMC are computationally intractable when you scale to hundreds or thousands of spatial units. Third, there's no existing bridge that brings vector surveillance data — the biting rates, mosquito survival, insecticide resistance data being collected by entomologists — into the routine case-data-driven mapping pipelines. Our framework addresses all three of these gaps.
-->
