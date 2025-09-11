
# ---- Packages ----
library(greta)
library(greta.dynamics)

# ---- Ross–Macdonald ODEs (vectorised-friendly for greta) ----
dx_dt <- function(x, z, m, a, b, r) {
  m * a * b * z * (1 - x) - r * x
}

dz_dt <- function(x, z, a, c, g) {
  a * c * x * (1 - z) - g * z
}

# ---- Constants as greta data (deterministic) ----
m  <- as_data(100)
a  <- as_data(0.1)
b  <- as_data(0.9)
r  <- as_data(1 / 21)   # human mean infectious duration = 21 days
c  <- as_data(0.9)
g  <- as_data(1 / 21)   # mosquito lifespan = 21 days
dt <- as_data(0.1)

x_init <- as_data(0.02)
z_init <- as_data(0.001)

# ---- Time grid & steps ----
t_max  <- 365
t      <- seq(0, t_max, by = as.numeric(dt))
n_step <- length(t) - 1  # number of Euler steps

# ---- Euler transition function for greta.dynamics ----
euler_step <- function(state, m, a, b, c, r, g, dt) {
  x <- state[1]
  z <- state[2]
  x_next <- x + dx_dt(x, z, m, a, b, r) * dt
  z_next <- z + dz_dt(x, z, a, c, g) * dt
  c(x_next, z_next)
}

# ---- Run the dynamics ----
solution <- iterate_dynamic_function(
  transition_function = euler_step,
  initial_state = c(x_init, z_init),
  niter = n_step,
  tol = 0,
  m = m, a = a, b = b, c = c, r = r, g = g, dt = dt
)

# Extract latent state time series (greta arrays)
x_ts <- solution$all_states[1, ]  # length n_step + 1 (includes initial state)
z_ts <- solution$all_states[2, ]

# ---- Evaluate to numeric (deterministic calculate) ----
sim <- calculate(x_ts, z_ts)

x_num <- as.numeric(sim$x_ts[1, 1, ])  # drop greta dims
z_num <- as.numeric(sim$z_ts[1, 1, ])

# ---- Plot (should match your base-R Euler loop) ----
op <- par(mfrow = c(1, 2))
plot(t, x_num, type = "l", xlab = "time (days)", ylab = "x (humans infected)")
plot(t, z_num, type = "l", xlab = "time (days)", ylab = "z (mosquitoes infected)")
par(op)
```

### Notes

# This is **deterministic**: we pass constants with `as_data()` and don’t put priors/likelihoods, so `calculate()` just evaluates the graph (no sampling).
# If you saw the earlier error about TensorFlow `2.16.0` vs `2.20.0`, align your TF runtime to **2.16.x** for current `greta` releases before running this.
# You can now layer a stochastic observation model on top (e.g., `poisson` likelihood on cases) without changing the dynamics block.
