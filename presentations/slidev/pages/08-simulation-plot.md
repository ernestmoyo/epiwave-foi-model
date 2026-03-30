---
layout: default
---

# Simulation Study: Rate Space vs Count Space

Two-panel view separating what the model predicts (rate) from what Stage 2 fits to (counts)

<div class="flex justify-center mt-2">
<img src="/images/sim_twopanel_site1.png" class="h-72 rounded shadow-lg" />
</div>

<div class="grid grid-cols-2 gap-4 mt-2 text-xs">

<div class="p-2 bg-blue-50 rounded border border-blue-200">

**Top — Rate Space:** Red dashed = I* (mechanistic rate from ODE). Blue = true incidence with GP residuals. The gap = exp(epsilon), showing how the GP distorts the mechanistic prediction. GP pushes true incidence above I* at some times, below at others.

</div>

<div class="p-2 bg-green-50 rounded border border-green-200">

**Bottom — Count Space:** Blue line = expected cases (gamma x I_true x N). Black dots = observed (Poisson draws). This is what Stage 2 actually fits to. The declining trend reflects ITN scale-up (0% to 70% coverage).

</div>

</div>

<!--
This two-panel view was designed to separate two different questions. The top panel shows rate space — the mechanistic prediction I-star from the ODE versus the true incidence including GP residuals. These are on the same scale, so you can see exactly where and how much the GP distorts the mechanistic prediction. The gap between the red dashed line and blue line IS the epsilon structure.

The bottom panel shows count space — what the model actually observes. Expected cases are gamma times I-true times population, and observed cases are Poisson draws around that expectation. Notice the cases are roughly 10% of the true incidence rate times population — that's the reporting rate gamma = 0.1.

Previously this was a single panel mixing rates and counts on the same axis, which made the mechanistic prediction invisible. The two-panel design was developed during our walkthrough session to show both quantities on appropriate scales.
-->
