# =========================================================
# Ross–Macdonald model in greta + greta.dynamics (end-to-end)
# =========================================================

# ---- Packages ----
library(greta)
library(greta.dynamics)

# ---------------------------------------------------------
# 1) Define the ODE right-hand sides (same equations as yours)
#    dx/dt and dz/dt are *deterministic* functions we’ll use
#    inside a simple Euler transition for greta.dynamics
# ---------------------------------------------------------
dx_dt <- function(x, z, m, a, b, r) m * a * b * z * (1 - x) - r * x
dz_dt <- function(x, z, a, c, g)   a * c * x * (1 - z) - g * z

# ---------------------------------------------------------
# 2) Constants as greta data (so they become greta arrays)
#    This mirrors your plain-R parameter choices
# ---------------------------------------------------------
m  <- as_data(100)      # mosquitoes per person
a  <- as_data(0.1)      # bites per mosquito per day
b  <- as_data(0.9)      # transmission: mosquito -> human
r  <- as_data(1 / 21)   # human recovery (1/duration)
c  <- as_data(0.9)      # transmission: human -> mosquito
g  <- as_data(1 / 21)   # mosquito loss rate (1/lifespan)
dt <- as_data(0.1)      # Euler step (days)

# Initial conditions (as greta arrays)
x_init <- as_data(0.02)   # fraction humans infected at t0
z_init <- as_data(0.001)  # fraction mosquitoes infected at t0

# ---------------------------------------------------------
# 3) Time grid (kept in base R for indexing/plotting)
# ---------------------------------------------------------
t_max  <- 365
t_seq  <- seq(0, t_max, by = as.numeric(dt))  # numeric time vector
n_step <- length(t_seq) - 1                   # number of Euler steps

# ---------------------------------------------------------
# 4) Define an Euler transition compatible with greta.dynamics
#    NOTE: signature must include `state` and `iter`.
#    Return the *next* state as a length-2 vector (x_next, z_next).
# ---------------------------------------------------------
euler_step <- function(state, m, a, b, c, r, g, dt, iter) {
  # unpack current state
  x <- state[1]
  z <- state[2]
  # one Euler step for each state component
  x_next <- x + dx_dt(x, z, m, a, b, r) * dt
  z_next <- z + dz_dt(x, z, a, c, g)   * dt
  c(x_next, z_next)
}

# ---------------------------------------------------------
# 5) Run the deterministic dynamics with greta.dynamics
#    iterate_dynamic_function builds a *greta graph* that
#    carries the full trajectory as greta arrays.
#    tol = 0 keeps it purely deterministic here.
# ---------------------------------------------------------
solution <- iterate_dynamic_function(
  transition_function = euler_step,
  initial_state = c(x_init, z_init),  # 2-vector (x0, z0)
  niter = n_step,                     # number of steps
  tol   = 0,                          # deterministic
  m = m, a = a, b = b, c = c, r = r, g = g, dt = dt
)

# `solution$all_states` is a 2 x (n_step + 1) greta array: [row 1 = x(t), row 2 = z(t)]
x_ts_greta <- solution$all_states[1, ]           # greta array (length n_step+1)
z_ts_greta <- solution$all_states[2, ]           # greta array (length n_step+1)

# ---------------------------------------------------------
# 6) Evaluate greta arrays to numeric with calculate()
#    calculate() returns a list; we name the elements for clarity.
# ---------------------------------------------------------
sim_vals <- calculate(list(x = x_ts_greta, z = z_ts_greta))
x_num <- as.numeric(sim_vals$x)   # numeric vector for plotting/use
z_num <- as.numeric(sim_vals$z)

# ---------------------------------------------------------
# 7) Quick checks analogous to your plain-R snippets
# ---------------------------------------------------------
# instantaneous rates at t0 (should match your dx_dt/dz_dt at init)
dx0 <- dx_dt(x = as.numeric(x_init), z = as.numeric(z_init), m = as.numeric(m),
             a = as.numeric(a), b = as.numeric(b), r = as.numeric(r))
dz0 <- dz_dt(x = as.numeric(x_init), z = as.numeric(z_init),
             a = as.numeric(a), c = as.numeric(c), g = as.numeric(g))
cat("dx/dt at t0:", dx0, "\n")
cat("dz/dt at t0:", dz0, "\n")

# one-step-ahead (dt = 0.1 day) expectation by formula vs greta trajectory
x_next_formula <- as.numeric(x_init) + dx0 * as.numeric(dt)
z_next_formula <- as.numeric(z_init) + dz0 * as.numeric(dt)
cat("x(t0 + dt) formula:", x_next_formula, " | greta:", x_num[2], "\n")
cat("z(t0 + dt) formula:", z_next_formula, " | greta:", z_num[2], "\n")

# ---------------------------------------------------------
# 8) Plot trajectories
# ---------------------------------------------------------
op <- par(mfrow = c(1, 2))
plot(t_seq, x_num, type = "l", xlab = "time (days)", ylab = "x (humans infected)",
     main = "Humans (x) via greta.dynamics")
plot(t_seq, z_num, type = "l", xlab = "time (days)", ylab = "z (mosquitoes infected)",
     main = "Mosquitoes (z) via greta.dynamics")
par(op)

