---
layout: default
---

# Stage 1: R Implementation

ODE solver with time-varying parameters and intervention effects

```r {1-8|10-16|18-21|all}{lines:true,maxHeight:'280px'}
ross_macdonald_ode <- function(t, state, parms) {
  x <- state[1]; z <- state[2]
  m <- if (is.function(parms$m)) parms$m(t) else parms$m
  a <- if (is.function(parms$a)) parms$a(t) else parms$a
  g <- if (is.function(parms$g)) parms$g(t) else parms$g
  b <- parms$b; c <- parms$c; r <- parms$r

  dx_dt <- m * a * b * z * (1 - x) - r * x
  dz_dt <- a * c * x * (1 - z) - g * z
  list(c(dx_dt, dz_dt))
}

# Mechanistic prediction: FOI x Population
compute_mechanistic_prediction <- function(m_matrix, a_matrix, b,
                                           z_matrix, population_matrix) {
  m_matrix * a_matrix * b * z_matrix * population_matrix
}

# Multi-site solver (one-time forward integration per site)
solve_ross_macdonald_multi_site <- function(m_matrix, a_matrix, g_matrix,
                                            times, b, c, r, ...) {
  for (site in 1:ncol(m_matrix)) {
    m_func <- approxfun(times, m_matrix[, site], rule = 2)
    # ... a_func, g_func similarly
    solution <- ode(y = c(x = 0.01, z = 0.001), times = times,
                    func = ross_macdonald_ode, parms = list(...))
  }
}
```

<div class="mt-1 text-xs opacity-70">

Uses `deSolve::ode()` with LSODA algorithm. Time-varying $m(t)$, $a(t)$, $g(t)$ via `approxfun()` interpolation. Interventions applied upstream via `apply_interventions()`.

</div>

<!--
Here's the actual R implementation, matching the cleaned codebase.

The ODE function takes state variables x and z, extracts time-varying parameters via function calls, then computes the Ross-Macdonald derivatives. Parameters m, a, g can be either constants or interpolation functions — this allows time-varying Vector Atlas data.

compute_mechanistic_prediction is a single line — it multiplies the Force of Infection components and scales by population to get I-star.

solve_ross_macdonald_multi_site loops over sites, creates approxfun interpolators for each site's parameter time series, and solves the ODE independently per site using LSODA.

The key point is that apply_interventions has already adjusted m, a, g before they reach the ODE — so ITN and IRS effects are baked into the time-varying parameters.
-->
