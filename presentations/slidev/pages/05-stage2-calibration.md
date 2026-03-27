---
layout: two-cols
---

# Stage 2: GP + Dual Likelihood

Gaussian Process residuals with mechanistic offset — **5 free parameters**

<div class="mt-2 text-sm">

**Latent infection incidence:**

$$\log(I_{s,t}) = \alpha + \log(I^*_{s,t}) + \varepsilon_{s,t}$$

**GP prior on residuals:**

$$\varepsilon \sim \text{GP}(0,\; \sigma^2 \cdot K_{\text{space}}(\phi) \cdot K_{\text{time}}(\theta))$$

**Dual likelihood:**

$$C_{s,t} \sim \text{Poisson}(\gamma \cdot I_{s,t})$$
$$Y_{s,t} \sim \text{Binomial}(N_{s,t},\; x_{s,t})$$

</div>

::right::

<div class="ml-4 mt-2">

**greta + greta.gp implementation:**

```r {lines:true}
fit_epiwave_gp <- function(
    observed_cases, I_star, x_star,
    coords, prev_data, ...) {

  # 5 free parameters
  alpha     <- normal(0, 1)
  gamma     <- normal(0.1, 0.05,
    truncation = c(0.001, Inf))
  log_sigma <- normal(-1, 1)
  log_phi   <- normal(1, 1)
  log_theta <- normal(1, 1)

  # GP residuals
  K <- build_gp_kernel(
    exp(log_sigma), exp(log_phi),
    exp(log_theta))
  epsilon <- gp(coords, K, ...)

  # Dual likelihood
  I <- exp(alpha + log(I_star) + epsilon)
  distribution(cases) <- poisson(gamma*I)
  distribution(prev) <- binomial(N, x_adj)
}
```

</div>

<!--
Stage 2 is where we calibrate to observed data. The structural equation uses log(I-star) — our mechanistic prediction from Stage 1 — as a fixed offset. The Gaussian Process epsilon captures spatially-correlated departures from the mechanistic model.

The key design choice is the dual likelihood. Poisson for case counts and Binomial for prevalence surveys. This is critical because with case counts alone, alpha and gamma are non-identifiable — they form a ridge in posterior space. But prevalence surveys are unaffected by case ascertainment, so they separately inform infection incidence. This makes both alpha and gamma individually identifiable.

The GP uses a separable kernel — Matern 5/2 for space times Exponential for time — following epiwave.mapping. With sparse GP inducing points, this scales to large spatial problems. The greta.gp package handles the kernel construction and GP sampling within the HMC framework.
-->
