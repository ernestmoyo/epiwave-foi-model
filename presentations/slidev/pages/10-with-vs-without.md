---
layout: center
---

# WITH vs WITHOUT Mechanistic Offset

<div class="flex justify-center">
<img src="/images/with_vs_without_offset.png" class="h-72 rounded shadow-lg" />
</div>

<div class="mt-2 text-xs text-center max-w-2xl mx-auto">

**What you're seeing:** The blue line (WITH offset) tracks the seasonal dynamics of observed cases (black dots) with RMSE = 14.7. The red dashed line (WITHOUT offset) collapses to the **grand mean** (~207 cases) — it has no temporal structure and RMSE = 125.5. The recovered reporting rate of 0.1003 (true = 0.10) confirms near-zero estimation bias. The mechanistic offset provides the temporal "shape" while $\lambda$ provides the "scale".

</div>

<!--
This is the key comparison plot. The blue curve — our WITH offset model — captures the full seasonal dynamics. It tracks the observed cases (black dots) through every peak and trough. RMSE is 14.7.

The red dashed line — the WITHOUT offset model — can only predict the grand mean. Without the mechanistic information telling it about seasonal dynamics, ITN effects, and spatial variation in vector parameters, it has no temporal structure. It just predicts roughly 207 cases every month. RMSE is 125.5.

The title bar shows the recovered reporting rate: 0.1003 versus the true 0.10. That's an 88% improvement in RMSE from including the mechanistic offset. This is the core result that demonstrates the value of incorporating vector data into malaria mapping.

Notice also that the blue line slightly overestimates at peaks and underestimates at troughs — this is because we're using the posterior mean of lambda, which applies a single scaling factor globally. In a real application with spatial random effects, this would be more flexible.
-->
