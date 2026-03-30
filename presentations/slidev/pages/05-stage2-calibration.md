---
layout: default
---

# Stage 2: GP + Dual Likelihood

Gaussian Process residuals with mechanistic offset — **5 free parameters**

<div class="grid grid-cols-2 gap-4 mt-2">

<div class="text-xs">

**Latent infection incidence:**

$$\log(I_{s,t}) = \alpha + \log(I^*_{s,t}) + \varepsilon_{s,t}$$

**Spatial GP + AR(1) temporal** (epiwave.mapping):

$$f_t \sim \text{GP}(0,\; \sigma^2 \cdot \text{Matérn}_{5/2}(\phi))$$
$$\varepsilon_t = \theta \cdot \varepsilon_{t-1} + f_t$$

**Dual likelihood:**

$$C_{s,t} \sim \text{Poisson}(\gamma \cdot I_{s,t} \cdot N_s)$$
$$Y_{s,t} \sim \text{Binomial}(N_{s,t},\; x_{s,t})$$

</div>

<div class="text-xs">

**What each parameter does:**

| Param | Prior | Role |
|:-----:|:------|:-----|
| $\alpha$ | N(0, 1) | Intercept — systematic bias in I\* |
| $\gamma$ | N(0.1, 0.05)+ | Reporting rate (case ascertainment) |
| $\sigma^2$ | LogN(-0.5, 0.5) | GP variance — how much reality departs from I\* |
| $\phi$ | LogN(0.5, 0.5) | Spatial lengthscale — range of spatial correlation |
| $\theta$ | Uniform(0, 1) | AR(1) temporal correlation |

<div class="mt-2 p-2 bg-green-50 rounded border border-green-200">

**Why dual likelihood?** With cases alone, $\alpha$ and $\gamma$ trade off (non-identifiable). Prevalence surveys inform $\alpha$ independently of $\gamma$.

</div>

<div class="mt-1 p-2 bg-blue-50 rounded border border-blue-200">

**Why GP + AR(1)?** Keeps latent space at n_sites (~10) not n_sites x n_times (~490). HMC works (<1% bad transitions).

</div>

</div>

</div>

<!--
Stage 2 is where we calibrate to observed data. The structural equation uses log(I-star) — our mechanistic prediction from Stage 1 — as a fixed offset. The Gaussian Process epsilon captures spatially-correlated departures from the mechanistic model.

The five parameters each have a clear role. Alpha adjusts the overall level — if alpha is 0, I-star is perfectly calibrated. Gamma is the reporting rate — what fraction of true infections become reported cases. Sigma-squared is the GP variance — how much spatial variation exists beyond what the mechanistic model explains. Phi is the spatial lengthscale — how far spatial correlation extends. Theta is the AR(1) correlation — how persistent residuals are across time.

The dual likelihood is critical. Poisson for case counts includes gamma (reporting rate), but Binomial prevalence surveys don't depend on gamma at all. So prevalence data separately informs alpha, making both parameters individually identifiable.

The spatial GP plus AR(1) pattern matches epiwave.mapping exactly. It keeps the latent space small — only n_sites variables per time step, not the full n_sites times n_times grid. This is what makes HMC tractable.
-->
