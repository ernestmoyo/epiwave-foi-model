---
layout: center
---

# Simulation Study — Stage 1 Output

<div class="flex justify-center">
<img src="/images/simulation_study_site1.png" class="h-72 rounded shadow-lg" />
</div>

<div class="mt-2 text-xs text-center max-w-2xl mx-auto">

**What you're seeing:** The blue line (true incidence) and red dashed line (mechanistic prediction $I^*$) overlap almost perfectly — confirming Stage 1 correctly recovers transmission dynamics. Black dots (observed cases) are ~10% of true incidence (reporting rate = 0.1), with Poisson noise. The seasonal oscillation reflects bimodal East African rainfall forcing on mosquito abundance, with declining amplitude from simulated ITN scale-up (0% → 70% coverage).

</div>

<!--
This plot shows the Stage 1 output for Site 1 in our simulation study. Three things to notice:

First, the blue line (true incidence) and the red dashed line (mechanistic prediction I-star) are essentially identical — they overlay each other. This confirms that when we solve the Ross-Macdonald ODEs with the correct fixed parameters, we recover the true transmission dynamics perfectly. This is expected in a simulation study where we know the true parameters.

Second, the black dots — observed cases — are much lower than the true incidence line. They represent about 10% of true incidence, which is our simulated reporting rate. There's also Poisson noise adding scatter.

Third, notice the seasonal pattern — two peaks per year from the bimodal East African rainfall seasonality, with declining amplitude over time. That declining trend is the simulated ITN scale-up from 0% to 70% coverage reducing transmission. The model correctly captures this intervention effect through the apply_interventions function modifying m, a, and g parameters.
-->
