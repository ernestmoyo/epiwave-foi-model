---
layout: two-cols
---

# Stage 2 Code: Dual Likelihood

<div class="text-xs mt-1">

```r
  # continued from fit_epiwave_gp()...

  # --- Poisson likelihood (cases) ---
  # Population enters HERE, not in I*
  pop_t <- t(population)
  cases_greta <- as_data(t(observed_cases))
  expected_cases <-
    gamma_rr * I_latent * pop_t
  distribution(cases_greta) <-
    poisson(expected_cases)

  # --- Binomial likelihood (prevalence) ---
  # Independent of gamma — breaks the ridge
  x_at_surveys <- x_star[survey_indices]
  log_odds_base <-
    log(x_at_surveys) - log(1 - x_at_surveys)
  log_odds_adj <-
    log_odds_base + epsilon[survey_indices]
  prev_prob <- ilogit(log_odds_adj)

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

Expected cases = gamma x I_latent x population

- **gamma** (reporting rate) scales the prediction
- **population** enters here, not in I\*
- This likelihood informs the product exp(alpha) x gamma, but cannot separate them alone

</div>

<div class="p-2 bg-green-50 rounded mb-3 border-l-3 border-green-400">

**Binomial likelihood (prevalence)**

Prevalence surveys measure actual infection in sampled populations. The GP residuals adjust the baseline ODE prevalence (x_star) on the logit scale.

- **Completely independent of gamma** — reporting rate does not affect survey results
- This is what separately informs alpha, breaking the alpha-gamma ridge
- Available at a subset of site-time combinations (30% in simulation)

</div>

<div class="p-2 bg-amber-50 rounded mb-3 border-l-3 border-amber-400">

**`model()` call**

Returns the 5 parameters for `greta::mcmc()` to sample via HMC. Everything else (GP, likelihoods, latent variables) is handled internally by TensorFlow's automatic differentiation.

</div>

</div>

<!--
This slide shows the dual likelihood — the key innovation in Stage 2. The Poisson likelihood for case counts includes population and gamma. The Binomial likelihood for prevalence surveys is independent of gamma. Together they make both alpha and gamma individually identifiable. The model call at the end returns all 5 parameters for MCMC sampling.
-->
