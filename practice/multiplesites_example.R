a <- c(1, 2, 3)
b <- pi

#c_i = a_i * b

c <- rep(NA, length(a))
for(i in seq_along(a)) {
  c[i] <- a[i] * b
}

c

c <- b * a
c



a <- c(1, 2, 3)
b <- c(1, 2)
c <- a * b

n_locations <- 3
n_times <- 10
m <- matrix(runif(n_locations * n_times),
            nrow = n_locations,
            ncol = n_times)
m

# dxi_dt = m_it * a * b * z_it * (1 - x_it) - r * x_it
for (i in seq_len(n_locations)) {
  for (t in seq_len(n_times)) {
    dxi_dt = m[i, t] * a * b * z[i, t] * (1 - x[i, t]) - r * x[i, t]
  }
}

# vectorised over locations
for (t in seq_len(n_times)) {
  dx_dt = m[, t] * a * b * z[, t] * (1 - x[, t]) - r * x[, t]
}

m[, 2]
