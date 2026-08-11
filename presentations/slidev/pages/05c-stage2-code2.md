---
layout: two-cols
---

# Stage 2 Code: Dual Likelihood <span class="text-xs align-middle px-2 py-0.5 rounded-full bg-emerald-100 text-emerald-700 font-bold">UPDATED · Aug 2026</span>

<div class="text-xs mt-1">

```r
  # continued from fit_epiwave_gp()...

  # --- Poisson likelihood (cases) ---
  # N_pop enters HERE, not in I*
  N_pop_t <- t(N_pop)
  cases_greta <- as_data(t(observed_cases))
  expected_cases <-
    gamma_rr * I_latent * N_pop_t
  distribution(cases_greta) <-
    poisson(expected_cases)

  # --- Binomial likelihood (prevalence) ---
  # Prevalence depends on I (GP-adjusted),
  # NOT mechanistic X -> breaks the ridge
  prev_intensity <-
    I_latent %*% prev_conv_matrix
  intensity <- prev_intensity[
    prev_data$survey_indices]
  prev_prob <- 1 - exp(-intensity)

  n_pos_greta <- as_data(n_positive)
  distribution(n_pos_greta) <-
    binomial(n_tested, prev_prob)

  # --- Model returns 5 parameters ---
  model(alpha, gamma_rr,
        sigma2, phi, theta)
}
```

</div>

::right::

<div class="ml-3 mt-6 text-xs">

<div class="p-2 bg-purple-50 rounded mb-3 border-l-3 border-purple-400">

**Poisson likelihood (cases)**

Expected cases = gamma x I_latent x **N_pop**

- **gamma** (reporting rate) scales the prediction
- **N_pop** (population) enters here, not in I\*

</div>

<div class="p-2 bg-green-50 rounded mb-3 border-l-3 border-green-400">

**Binomial likelihood (prevalence)**

Prevalence now derived from **I** (the GP-adjusted incidence) via the detectability convolution: `intensity = I convolved with q`, then `prev = 1 - exp(-intensity)`.

- Because I carries **alpha + epsilon**, this likelihood **depends on alpha** — separately informing it and breaking the ridge.

</div>

<div class="p-2 bg-amber-50 rounded border-l-3 border-amber-400">

**What changed since March — the key fix**

In March, prevalence was computed from x_star, the mechanistic prevalence, on the logit scale. Because x_star does not contain alpha, it could not separate alpha from gamma. It now depends on the modelled incidence I through the detectability convolution (Nick).

</div>

</div>

<!--
This is the slide with the most important code change. The Poisson likelihood is the same idea — expected cases equal gamma times latent incidence times population, now called N_pop. The big change is the prevalence likelihood. In March it was written against x_star, the mechanistic prevalence from the ODE, on the logit scale. The problem, which Nick identified, is that x_star does not depend on alpha, so that likelihood could never break the alpha–gamma ridge. The current code derives prevalence from I, the GP-adjusted latent incidence, by convolving it with the test-detectability kernel and mapping through one minus exp minus intensity. Because I carries alpha and the GP residual, the prevalence data now depends on alpha and separately informs it. That is the fix behind the 80 percent coverage result.
-->
