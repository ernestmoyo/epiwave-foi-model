---
layout: two-cols
---

# Stage 2 Code: Priors and GP

<div class="text-xs mt-1">

```r
fit_epiwave_gp <- function(
    observed_cases, I_star, x_star,
    population, spatial_coords,
    prev_data, use_mechanistic = TRUE) {

  # --- 5 free parameters ---
  alpha    <- normal(0, 1)
  gamma_rr <- normal(0.1, 0.05,
                truncation = c(0.001, Inf))
  sigma2   <- lognormal(-0.5, 0.5)
  phi      <- lognormal(0.5, 0.5)
  theta    <- variable(lower=0, upper=1)

  # --- Spatial GP + AR(1) temporal ---
  K <- build_gp_kernel(phi, sigma2)
  f <- gp(spatial_coords, K,
           n = n_times, tol = 1e-3)
  epsilon_mat <- ar1(rho = theta,
                     innovations = f)

  # --- Latent incidence ---
  if (use_mechanistic) {
    log_I <- alpha + log(I_star)
             + epsilon_mat
  } else {
    log_I <- alpha + epsilon_mat
  }
  I_latent <- exp(log_I)
```

</div>

::right::

<div class="ml-3 mt-6 text-xs">

<div class="p-2 bg-blue-50 rounded mb-3 border-l-3 border-blue-400">

**5 greta parameters**

All defined as greta variables — TensorFlow computes gradients for HMC automatically. `alpha` (intercept), `gamma_rr` (reporting rate), `sigma2` (GP variance), `phi` (spatial lengthscale), `theta` (AR1 correlation).

</div>

<div class="p-2 bg-green-50 rounded mb-3 border-l-3 border-green-400">

**Spatial GP + AR(1)**

`gp()` draws spatial innovations per time step using Matern 5/2. `ar1()` correlates them across time. Latent space = **n_sites only** (~10 variables), not n_sites x n_times (~490). This is what makes HMC tractable.

</div>

<div class="p-2 bg-gray-100 rounded border-l-3 border-gray-400">

**`use_mechanistic` toggle**

TRUE: log(I) = α + log(I\*) + ε (WITH offset)
FALSE: log(I) = α + ε (standard geostatistical)

</div>

</div>

<!--
This slide shows the first half of the Stage 2 fitting function — the parameter definitions and GP structure. Five greta variables with their priors, the spatial GP plus AR(1) temporal correlation, and the latent incidence calculation. The use_mechanistic flag toggles the core comparison.
-->
