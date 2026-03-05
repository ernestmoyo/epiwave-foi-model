---
layout: two-cols
---

# Stage 2: Bayesian Calibration

Negative Binomial likelihood with mechanistic offset — **only 2 free parameters**

<div class="mt-2 text-sm">

**Structural Equation:**

$$\log(\mu_{s,t}) = \lambda + \log(I^*_{s,t})$$

**Likelihood:**

$$C_{s,t} \sim \text{NegBin}(\text{size},\; p_{s,t})$$

where $p_{s,t} = \frac{\text{size}}{\text{size} + \mu_{s,t}}$

**Priors:**

$$\lambda \sim \text{Normal}(-2, 1)$$
$$\log(\text{size}) \sim \text{Normal}(3, 1)$$
$$\text{size} = \exp(\log\text{\_size})$$

</div>

::right::

<div class="ml-4 mt-2">

**greta implementation:**

```r {lines:true}
fit_epiwave_with_offset <- function(
    observed_cases, I_star,
    use_mechanistic = TRUE) {

  library(greta)
  cases_vec  <- as.vector(observed_cases)
  I_star_vec <- as.vector(I_star)

  # Priors — only 2 free parameters
  log_rate <- normal(-2, 1)
  log_size <- normal(3, 1)
  size     <- exp(log_size)

  # Linear predictor with offset
  log_mu <- log_rate + log(I_star_vec + 1e-10)
  mu     <- exp(log_mu)
  prob   <- size / (size + mu)

  # Likelihood
  distribution(cases_vec) <-
    negative_binomial(size, prob)

  model(log_rate, log_size)
}
```

</div>

<!--
Stage 2 is where we calibrate to observed case data. The structural equation uses log(I-star) — our mechanistic prediction from Stage 1 — as a fixed offset. The parameter lambda captures the log reporting rate, which is the only identifiable scaling parameter.

We use a Negative Binomial likelihood because it naturally handles overdispersion. The log_size parameter controls overdispersion on the log scale — when size is large, we recover the Poisson as a special case. We parameterise on the log scale with Normal(3,1) to give HMC a smooth, unconstrained posterior surface.

The greta implementation on the right shows how compact this is — the entire model specification is just a few lines. We define two priors, construct the linear predictor with the mechanistic offset, and specify the likelihood. That's it — ready for HMC sampling.

The Normal(3,1) prior on log_size gives a median size of ~20, with 95% CI from ~1 to ~400 — broad but regularising. This replaced an earlier Beta(1,9) prior on phi that caused boundary issues with ESS=82.
-->
