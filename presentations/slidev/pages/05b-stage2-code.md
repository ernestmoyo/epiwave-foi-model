---
layout: two-cols
---

# Stage 2 Code: Priors and GP <span class="text-xs align-middle px-2 py-0.5 rounded-full bg-emerald-100 text-emerald-700 font-bold">UPDATED · Aug 2026</span>

<div class="text-xs mt-1">

```r
fit_epiwave_gp <- function(
    observed_cases, I_star, N_pop,
    spatial_coords, prev_data = NULL,
    prev_conv_matrix = NULL,
    use_mechanistic = TRUE,
    inducing = NULL, gp_tol = 1e-3,
    center_alpha = FALSE) {

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
           n = n_times, tol = gp_tol)
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

`alpha` (intercept), `gamma_rr` (reporting rate), `sigma2` (GP variance), `phi` (spatial lengthscale), `theta` (AR1). TensorFlow computes HMC gradients automatically.

</div>

<div class="p-2 bg-green-50 rounded mb-3 border-l-3 border-green-400">

**Spatial GP + AR(1)**

`gp()` draws spatial innovations per time step (Matern 5/2); `ar1()` correlates across time. Latent space = **n_sites only** (~10), not n_sites x n_times — keeps HMC tractable.

</div>

<div class="p-2 bg-amber-50 rounded border-l-3 border-amber-400">

**What changed since March**

Signature now takes **`N_pop`** (was `population`) and a **`prev_conv_matrix`**; the old **`x_star`** argument is gone — prevalence no longer uses mechanistic X (see next slide).

</div>

</div>

<!--
This is the first half of the Stage 2 fitting function, updated to the current code. The parameter block and GP structure are unchanged — five greta variables, spatial GP plus AR(1). Two things in the signature changed since March: population is now called N_pop, to keep it distinct from the survey sample size, and there is a new prev_conv_matrix argument. The old x_star argument is gone, because prevalence no longer depends on the mechanistic prevalence X — that is the fix on the next slide.
-->
