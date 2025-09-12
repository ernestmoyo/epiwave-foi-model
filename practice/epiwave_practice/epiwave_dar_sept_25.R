## ------------------------------------------------------------
## EpiWave-FOI: minimal per-pixel simulator (discrete time)
## Two state variables:
##   x_t = fraction of humans infected
##   z_t = fraction of mosquitoes infected
##
## Flows (per day):
##   lambda_t = m_t * a_t * b * z_t             (human FOI)
##   incidence_t = lambda_t * (1 - x_t)         (new human infections per capita)
##   dx/dt ≈ lambda_t*(1 - x_t) - r*x_t
##   dz/dt ≈ a_t*c*x_t*(1 - z_t) - g*z_t
##
## We update with Euler steps:
##   x_{t+1} = x_t + dt * [ lambda_t*(1 - x_t) - r*x_t ]
##   z_{t+1} = z_t + dt * [ a_t*c*x_t*(1 - z_t) - g*z_t ]
##
## All rates are per day; states are proportions in [0, 1].
## ------------------------------------------------------------

simulate_epi <- function(
    days = 365,     # number of days to simulate
    dt   = 1,       # day step (keep at 1 for clarity)
    # ---- fixed biological/epi parameters (tune as needed) ----
    b = 0.3,        # Pr(mosq -> human) per bite
    c = 0.5,        # Pr(human -> mosq) per bite
    r = 1/200,      # human recovery rate (mean infectious ~200 days)
    g = 1/10,       # mosquito death rate (mean lifespan ~10 days)
    # ---- biting & abundance (can be constant or seasonal) ----
    m0 = 5.0,       # baseline mosquitoes per human
    a0 = 0.3,       # baseline bites per mosquito per day
    season_m_amp = 0.3,   # amplitude for m seasonality (0 = constant)
    season_a_amp = 0.1,   # amplitude for a seasonality (0 = constant)
    season_phase = 0,     # phase shift in radians (align peaks)
    # ---- initial conditions (fractions) ----
    x0 = 0.01,      # 1% humans infected initially
    z0 = 0.02       # 2% mosquitoes infected initially
) {
  # Pre-allocate output vectors for speed/clarity
  T <- days + 1L                  # include t = 0
  time <- seq(0, days, by = 1)
  x <- numeric(T); z <- numeric(T)
  lambda <- numeric(T)
  inc <- numeric(T)               # per-capita incidence rate (per day)
  m_t <- numeric(T); a_t <- numeric(T)

  # Helper: daily seasonality using a sine wave on [0, days]
  #   value_t = base * (1 + amp * sin(2*pi*t/365 + phase)), constrained positive
  seasonal <- function(base, amp, t, phase = 0) {
    val <- base * (1 + amp * sin(2*pi*t/365 + phase))
    pmax(val, 1e-9) # ensure strictly positive
  }

  # Set initial states
  x[1] <- x0
  z[1] <- z0

  # Main simulation loop (Euler updates)
  for (i in 1:T) {
    t <- time[i]

    # Compute time-varying m_t and a_t (seasonal or constant)
    m_t[i] <- seasonal(m0, season_m_amp, t, season_phase)
    a_t[i] <- seasonal(a0, season_a_amp, t, season_phase)

    # Current FOI on humans (per susceptible per day)
    lambda[i] <- m_t[i] * a_t[i] * b * z[i]

    # Per-capita incidence (new human infections per person per day)
    inc[i] <- lambda[i] * (1 - x[i])

    # Skip last step’s update
    if (i == T) break

    # Euler updates for x and z
    dx <- dt * ( inc[i] - r * x[i] )
    dz <- dt * ( a_t[i] * c * x[i] * (1 - z[i]) - g * z[i] )

    x[i + 1] <- x[i] + dx
    z[i + 1] <- z[i] + dz

    # Keep states in [0, 1] to avoid numerical drift
    x[i + 1] <- min(max(x[i + 1], 0), 1)
    z[i + 1] <- min(max(z[i + 1], 0), 1)
  }

  # Bundle results
  out <- data.frame(
    day = time,
    x = x,                     # human infected fraction
    z = z,                     # mosquito infected fraction
    lambda = lambda,           # FOI on humans (per day)
    incidence_per_capita = inc,# new infections per person per day
    m = m_t,                   # mosquitoes per human (time-varying)
    a = a_t                    # bites per mosquito per day (time-varying)
  )
  return(out)
}

## ------------------------------------------------------------
## Example run
## ------------------------------------------------------------

sim <- simulate_epi(
  days = 365,
  # Try zero seasonality first for a clean equilibrium feel:
  season_m_amp = 0.0,
  season_a_amp = 0.0
)

head(sim, 10)   # peek at the first 10 rows

## ------------------------------------------------------------
## Simple base-R plots (optional)
## ------------------------------------------------------------

# Humans infected (x) and mosquitoes infected (z)
plot(sim$day, sim$x, type = "l", xlab = "Day", ylab = "Fraction",
     main = "States: Humans (x) and Mosquitoes (z)")
lines(sim$day, sim$z, lty = 2)
legend("topright", legend = c("x (humans infected)", "z (mosquitoes infected)"),
       lty = c(1,2), bty = "n")

# FOI and incidence (per capita)
plot(sim$day, sim$lambda, type = "l", xlab = "Day", ylab = "Rate per day",
     main = "FOI (lambda) and Per-capita Incidence")
lines(sim$day, sim$incidence_per_capita, lty = 2)
legend("topright", legend = c("lambda (FOI)", "incidence per capita"),
       lty = c(1,2), bty = "n")

# If you want to see seasonality effects, re-run with non-zero amplitudes:
# sim_season <- simulate_epi(season_m_amp = 0.3, season_a_amp = 0.1)
# plot(sim_season$day, sim_season$incidence_per_capita, type="l",
#      xlab="Day", ylab="Rate per day",
#      main="Incidence with seasonal m and a")
