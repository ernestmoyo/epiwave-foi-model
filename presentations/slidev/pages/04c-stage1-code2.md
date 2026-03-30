---
layout: two-cols
---

# Stage 1 Code: I* and Interventions

<div class="text-xs mt-1">

```r
# I* is a RATE — no population here
compute_mechanistic_prediction <- function(
    m_matrix, a_matrix, b, z_matrix) {
  m_matrix * a_matrix * b * z_matrix
}

# ITN/IRS effects on entomological params
apply_interventions <- function(
    m, a, g,
    itn_coverage = NULL,
    irs_coverage = NULL,
    resistance_index = 0) {
  u <- 1 - resistance_index
  m_adj <- m; a_adj <- a; g_adj <- g

  if (!is.null(itn_coverage)) {
    # Reduce mosquito abundance
    m_adj <- m_adj * ((1 - itn_coverage)
      + itn_coverage*(1 - u*kill_rate))
    # Reduce biting rate
    a_adj <- a_adj *
      (1 - itn_coverage*u*feeding_inhibit)
    # Increase mortality
    g_adj <- g_adj *
      (1 + itn_coverage*u*mortality_boost)
  }
  # IRS similarly adjusts g and a

  list(m = m_adj, a = a_adj, g = g_adj)
}
```

</div>

::right::

<div class="ml-3 mt-6 text-xs">

<div class="p-2 bg-amber-50 rounded mb-3 border-l-3 border-amber-400">

**`compute_mechanistic_prediction()`**

One line: **I\* = m x a x b x z**

I\* is a **rate**, not a count. Population enters the Poisson likelihood in Stage 2, not here. This is the fixed offset that Stage 2 calibrates against.

</div>

<div class="p-2 bg-purple-50 rounded mb-3 border-l-3 border-purple-400">

**`apply_interventions()`**

Adjusts m, a, g **before** the ODE runs:

- **ITN (bed nets):** reduces mosquito abundance (m), reduces biting rate (a), increases mortality (g)
- **IRS (indoor spraying):** increases mortality (g), slightly reduces biting (a)
- **Resistance index** (0–1) scales intervention effectiveness down

This means intervention counterfactuals only require re-running Stage 1 — Stage 2 stays the same.

</div>

</div>

<!--
The mechanistic prediction is a single line — m times a times b times z. I-star is a rate — population enters the Poisson likelihood in Stage 2.

The apply_interventions function is important for the operational value of the framework. It adjusts entomological parameters upstream of the ODE based on ITN and IRS coverage. This means if you want to model a counterfactual — what happens if ITN coverage increases from 30% to 70% — you just re-run Stage 1 with new coverage values. Stage 2 stays the same. The resistance index allows modelling the impact of insecticide resistance on intervention effectiveness.
-->
