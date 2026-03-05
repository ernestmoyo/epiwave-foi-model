---
layout: default
---

# MCMC Convergence Diagnostics

<div class="flex justify-center">
<img src="/images/mcmc_trace_plots.png" class="h-40 rounded shadow-lg" />
</div>

<div class="grid grid-cols-2 gap-2 mt-1">

<div class="p-2 rounded-lg bg-blue-50 border border-blue-200 text-xs">

**WITH Mechanistic Offset**

| Parameter | $\hat{R}$ | ESS | Status |
|:---------:|:----------:|:---:|:------:|
| log_rate ($\lambda$) | 1.004 | 1,095 | Converged |
| log_size | 1.004 | 489 | Converged |

<div class="mt-1 text-xs opacity-80">
log_size: mean=5.49, sd=0.49 (size≈240). Reparam from Beta(phi) → Normal(log_size) — ESS 82→<b>489</b>.
</div>

</div>

<div class="p-2 rounded-lg bg-amber-50 border border-amber-200 text-xs">

**WITHOUT Mechanistic Offset**

| Parameter | $\hat{R}$ | ESS | Status |
|:---------:|:----------:|:---:|:------:|
| log_rate ($\lambda$) | 1.002 | 1,540 | Converged |
| phi ($\phi$) | 1.000 | 2,769 | Converged |

<div class="mt-1 text-xs opacity-80">
phi=0.357 (size≈1.8) — heavy overdispersion, no mechanistic structure absorbs variation.
</div>

</div>

</div>

<div class="mt-1 p-1 bg-gray-100 text-xs rounded text-center">

**Overdispersion Story:** WITH size≈240 (near-Poisson) vs WITHOUT size≈1.8 (heavy overdispersion) — mechanistic model explains virtually all variation

</div>

<!--
These are the MCMC trace plots for both models. The top row shows the WITH offset model, the bottom shows WITHOUT.

For log_rate, both models converge well — R-hat near 1, good ESS. The WITH model estimates log_rate around minus 2.30, which translates to a reporting rate of about 10% — matching our true value.

The interesting story is the overdispersion parameter. After reparameterising from phi ~ Beta(1,9) to log_size ~ Normal(3,1), the WITH offset model now converges cleanly — ESS jumped from 82 to 489. The estimated log_size of 5.49 translates to size ~240, which is near-Poisson. This confirms the mechanistic offset already explains the temporal and spatial variation.

In the WITHOUT offset model, phi is 0.357, meaning size is about 1.8 — very heavy overdispersion. Without the mechanistic structure, all that temporal variation looks like overdispersion to the model.

The reparameterisation lesson: when a parameter gets pushed to a boundary by the data, move to an unconstrained scale. The log transform solved the boundary mixing problem completely.
-->
