---
layout: default
---

# Stage 1: R Implementation

ODE solver with time-varying parameters and intervention effects

```r {1-5|7-14|16-19|all}{lines:true,maxHeight:'300px'}
ross_macdonald_ode <- function(t, state, parms) {
  x <- state[1]  # Human infection prevalence
  z <- state[2]  # Mosquito infection prevalence
  m <- if (is.function(parms$m)) parms$m(t) else parms$m
  a <- if (is.function(parms$a)) parms$a(t) else parms$a
  g <- if (is.function(parms$g)) parms$g(t) else parms$g
  b <- parms$b; c <- parms$c; r <- parms$r

  # Ross-MacDonald differential equations
  dx_dt <- m * a * b * z * (1 - x) - r * x
  dz_dt <- a * c * x * (1 - z) - g * z
  return(list(c(dx_dt, dz_dt)))
}

# Mechanistic prediction: FOI × Population
compute_mechanistic_prediction <- function(m, a, b, z, population) {
  foi <- m * a * b * z           # Force of Infection per person per day
  I_star <- foi * population     # Expected incidence
  return(I_star)
}
```

<div class="mt-1 text-xs opacity-70">

Uses `deSolve::ode()` with LSODA algorithm — automatically switches between stiff/non-stiff methods. Time-varying $m(t)$, $a(t)$, $g(t)$ via `approxfun()` interpolation.

</div>

<!--
Here's the actual R implementation. Click through to see the code build up.

First, the ODE function itself — it takes the state variables x and z, extracts parameters which can be functions of time for time-varying entomological parameters, then computes the derivatives.

The parameters m, a, g can be either constants or interpolation functions — this allows us to feed in time-varying Vector Atlas data.

Then the mechanistic prediction simply multiplies the Force of Infection components together and scales by population.

We use the LSODA algorithm from deSolve which automatically handles stiff and non-stiff regions of the ODE system. The time-varying parameters are handled via approxfun interpolation between data points.
-->
