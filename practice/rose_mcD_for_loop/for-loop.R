# Ross-Macdonald Differential Equations
dx_dt <- function(x, z, m, a, b, r) {
  m * a * b * z * (1 - x) - r * x
}

dz_dt <- function(x, z, a, c, g) {
  a * c * x * (1 - z) - g * z
}

# dt = 1 day

m <- 100
a <- 0.1
b <- 0.9
infection_duration <- 21
r <- 1/infection_duration
c <- 0.9
mosquito_lifespan <- 21
g <- 1 / mosquito_lifespan


x_init <- 0.02
z_init <- 0.001

# change in fraction of people infected per unit time (one day) from initial
# conditions
dx_dt(x = x_init, z = z_init, m = m, a = a, b = b, r = r)

# change in fraction of mosquitoes infected per unit time (one day) from initial
# conditions
dz_dt(x = x_init, z = z_init, a = a, c = c, g = g)

# fraction of people infected 1 day in the future
x_next_day <- x_init + dx_dt(x = x_init, z = z_init, m = m, a = a, b = b, r = r)

# fraction of mosquitoes infected 1 day in the future
z_next_day <- z_init + dz_dt(x = x_init, z = z_init, a = a, c = c, g = g)


# fraction of people infected 2 days in the future
x_in_two_days <- x_init + 2 * dx_dt(x = x_init, z = z_init, m = m, a = a, b = b, r = r)

# fraction of people infected 0.1 days in the future
x_in_0.1_days <- x_init + 0.1 * dx_dt(x = x_init, z = z_init, m = m, a = a, b = b, r = r)


t_max <- 10
dt <- 0.1
t <- seq(0, t_max, by = dt)
n_t <- length(t)
m <- 8 + (100 - 8)*(t/t_max) # in general this has to be close to a sine graph (this illustrative ie a linear ramp)

x <- rep(NA, n_t)
z <- rep(NA, n_t)
x[1] <- x_init
z[1] <- z_init
for(i in 2:n_t) {
  x[i] <- x[i - 1] + dx_dt(x = x[i - 1],
                           z = z[i - 1],
                           m = m[i - 1],
                           a = a,
                           b = b,
                           r = r) * dt
  z[i] <- z[i - 1] + dz_dt(x = x[i - 1],
                           z = z[i - 1],
                           a = a,
                           c = c,
                           g = g) * dt
}
par(mfrow = c(1, 2))
plot(x ~ t, type = "l")
plot(z ~ t, type = "l")
