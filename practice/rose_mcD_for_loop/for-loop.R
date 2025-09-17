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


# Multisite extension
s <- 10000 # number of sites - which can be changed as needed


t_max <- 10
dt <- 0.1
t <- seq(0, t_max, by = dt)
n_t <- length(t)
# m <- 8 + (100 - 8)*(t/t_max) # in general this has to be close to a sine graph (this illustrative ie a linear ramp)

# Initial conditions per site
x0 <- rep(x_init, s) # E.G. , c(0.02, 0.015, 0.03)
z0 <- rep(z_init, s)

# Site specific, time varying m(t) ILLUSTRATIVE LINEAR RAMP
set.seed(20250916)
m0 <- runif(s, 0, 50) # start value per site
mT <- runif(s, 50, 100) # end value per site
# m_mat: S x n_t (ROWS = SITES, COL = TIME)
m_mat <- sapply(t, function(tt) m0 + (mT - m0)*(tt/t_max))

# state matrices (rows = sites, cols = time)
x_mat <- matrix(NA_real_, nrow = s, ncol = n_t)
z_mat <- matrix(NA_real_, nrow = s, ncol = n_t)
x_mat[, 1] <- x0
z_mat[, 1] <- z0

# Single time loop, vector operations accross sites
for (j in 2:n_t) {
  x_prev <- x_mat[, j-1]
  z_prev <- z_mat[, j-1]
  m_now <- m_mat[, j-1]

  # a, b, r, c, g are shared across sites - making them length s vectors if needed
  dx <- dx_dt(x = x_prev, z = z_prev, m = m_now, a = a, b = b, r = r)
  dz <- dz_dt(x = x_prev, z = z_prev,            a = a, c = c, g = g)

  x_mat[, j] <- x_prev + dx * dt
  z_mat[ ,j] <- z_prev + dz * dt
}

# plotting one line per sire
par(mfrow = c(1, 2))
matplot(t, t(x_mat), type = "l")
matplot(t, t(z_mat), type = "l")
