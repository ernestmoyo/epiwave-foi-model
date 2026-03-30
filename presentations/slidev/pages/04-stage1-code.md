---
layout: two-cols
---

# Stage 1 Code: ODE System

<div class="text-xs mt-1">

```r
ross_macdonald_ode <- function(t, state, parms) {
  x <- state[1]; z <- state[2]
  # m, a, g: constant OR time-varying function
  m <- if (is.function(parms$m)) parms$m(t)
       else parms$m
  a <- if (is.function(parms$a)) parms$a(t)
       else parms$a
  g <- if (is.function(parms$g)) parms$g(t)
       else parms$g
  b <- parms$b; c <- parms$c; r <- parms$r
  dx_dt <- m*a*b*z*(1-x) - r*x
  dz_dt <- a*c*x*(1-z) - g*z
  list(c(dx_dt, dz_dt))
}

solve_ross_macdonald_multi_site <- function(
    m_matrix, a_matrix, g_matrix,
    times, b, c, r) {
  for (site in 1:ncol(m_matrix)) {
    m_func <- approxfun(times,
      m_matrix[, site], rule = 2)
    # a_func, g_func similarly
    solution <- ode(
      y = c(x = 0.01, z = 0.001),
      times = times,
      func = ross_macdonald_ode,
      parms = list(...),
      method = "lsoda")
  }
}
```

</div>

::right::

<div class="ml-3 mt-6 text-xs">

<div class="p-2 bg-blue-50 rounded mb-3 border-l-3 border-blue-400">

**`ross_macdonald_ode()`**

The ODE system. Parameters m, a, g accept either constants or `approxfun()` interpolators — so time-varying Vector Atlas data plugs in directly. Returns derivatives dx/dt, dz/dt.

</div>

<div class="p-2 bg-green-50 rounded mb-3 border-l-3 border-green-400">

**`solve_ross_macdonald_multi_site()`**

Loops over sites independently. Each site gets `approxfun()` interpolators for its m(t), a(t), g(t) time series, then solves the ODE **once** using LSODA.

Not called at every MCMC iteration — solved once up front. `apply_interventions()` adjusts m, a, g **before** this function runs.

</div>

</div>

<!--
This slide shows the first two Stage 1 functions. The ODE function handles time-varying parameters — if the parameter is a function from approxfun, it evaluates at time t; otherwise it uses a constant. The multi-site solver loops over spatial locations independently, solving each once up front.
-->
