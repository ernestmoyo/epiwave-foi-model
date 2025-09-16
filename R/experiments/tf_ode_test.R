# Test out different options for ODE solvers in tensorflow

# Steps:

# 1. Set up a moderately-realistic version of the Ross-MacDonald ODEs with
# time-varying parameters (strong bimodal seasonality, over 4 years) and N
# sites.

# 2. Write Tensorflow code for spline interpolation of monthly timeseries data
# to evaluate as a continuous function.

# 3. Write Tensorflow code to solve the ODEs with some parameters as continuous
# functions of time.

# 4. Speed-test the forward-mode solution of these ODEs with different solvers,
# for increasing number of sites N.

# 5. Speed-test differentiation of these ODEs with different solvers,
# for increasing number of sites N.

# 6. Wrap the most optimal solver up as a greta operation, for use in a greta
# model.



# load packages

# load greta first (even though we aren't using it) so that we are using the
# same version of tensorflow
library(greta)
library(tensorflow)
library(deSolve)


# 1. Set up a moderately-realistic version of the Ross-MacDonald ODEs with
# time-varying parameters (strong bimodal seasonality, over 4 years) and N
# sites.

# These are the base equations:

# dx/dt = m a b z (1 - x) - r x
# dz/dt = a c x (1 - z) - g z

# We want to make some of these these time-varying, so we will define continuous
# functions of time m(t), a(t), and g(t), since we will
# model these using entomological and vector control data:

# dx/dt = m(t) a(t) b z (1 - x) - r x
# dz/dt = a(t) c x (1 - z) - g(t) z

# We want to solve this simultaneously across N sites, each site i with
# different functions: m_i(t), a_i(t) and g_i(t)), and for B 'batches' b
# (different simulation runs, or MCMC chains), with independent timeseries
# within each site and batch. So we have:

# dx_{b,i}/dt = m_i(t) a_i(t) b z_{b,i} (1 - x_{b,i}) - r x_{b,i}
# dz_{b,i}/dt = a_i(t) c x_{b,i} (1 - z_{b,i}) - g_i(t) z

# Though we will compute these in a vectorised manner.

# We want to have seasonality in the time-varying parameters. This will lead to
# changes in prevalence over time, which may frustrate the ODE solvers. We need
# the test case to have sufficiently strong intra-annual seasonality to capture
# the likely changes in a real analysis. The most extreme cases (rapid changes
# in prevalence) will either be moderate but bimodal seasonality (as in e.g.
# Uganda, Tanzania), or strong unimodal seasonality (e.g. as in the Sahel). So
# we will combine these, to have strong bimodal seasonality (slightly
# unrealistically hard):

# make a curve with bimodal intra-annual seasonality, one peak higher than other
intra_annual_seasonality <- function(day) {
  year <- day / 365
  cycle <- 2 * pi
  sin(2 * cycle * year) - 0.5 * cos(cycle * year)
}

# and another with inter-annual seasonality, going up for 5 years then down for
# 5, with half as much influence as the annual variation
inter_annual_seasonality <- function(day) {
  year <- day / 365
  decade <- year / 10
  cycle <- 2 * pi
  -0.1 * cos(cycle * decade)
}

# combine these two
seasonality <- function(day) {
  intra_annual_seasonality(day) +
    inter_annual_seasonality(day)
}

# define m a and g functions

# m must be positive, with large integers
m <- function(day) {
  s <- seasonality(day)
  exp(1 + 0.5 * s)
}

# we are assuming biting rate is between 0 and 1, make biting activity slightly
# higher in the peak season (higher temperature and humidity)
a <- function(day) {
  s <- seasonality(day)
  plogis(-2 + 0.5 * s)
}

g <- function(day) {
  # model average mosquito lifespan based on seasonality
  s <- seasonality(day)
  lifespan <- exp(2 + 0.5 * s)
  # instantaeous death rate is 1 divided by lifespan
  1 / lifespan
}

# plot these to check

n_years <- 3
day <- seq(0, n_years * 365, length.out = 1e3)

par(mfrow = c(2, 2))
plot(seasonality(day) ~ day,
     type = "l",
     main = "seasonality")
plot(m(day) ~ day,
     type = "l",
     main = "abundance")
plot(a(day) ~ day,
     type = "l",
     main = "biting rate",
     ylim = c(0, 1))
plot(g(day) ~ day,
     type = "l",
     main = "mortality")

# solve these ODEs in deSolve to check we can do this correctly
deriv <- function(t, state, params) {
  with(as.list(c(state, params)), {
    m_t <- m(t)
    a_t <- a(t)
    g_t <- g(t)

    dx_dt <- m_t * a_t * b * z * (1 - x) - r * x
    dz_dt <- a_t * c * x * (1 - z) - g_t * z

    return(list(c(dx_dt, dz_dt)))
  })
}

# fix the other parameters, for our ODE solving experiments
params <- list(
  b = 0.8,
  c = 0.8,
  r = 1 / 7
)

# burn in for 1 year prior to simulations (start at 0 months, evaluate from 12
# months onwards)
days_burnin <- 365 * 5
days_sim <- 365 * 4
days_total <- days_burnin + days_sim
days <- seq(0, days_total, by = 1)
init <- c(x = 0.1, z = 0.01)
out <- ode(init, days, deriv, params)

# crop out the burnin, and relabel the months
keep <- days > days_burnin
out_sim <- out[keep, ]
out_sim[, "time"] <- out_sim[, "time"] - days_burnin

# compute the force of infection
foi <- m(out_sim[, "time"]) *
  a(out_sim[, "time"]) *
  params$b *
  out_sim[, "z"]

# and the infection incidence (per month, as a fraction of the population)
inf_inc <- foi * (1 - out_sim[, "x"])

par(mfrow = c(2, 2))
plot(x ~ time,
     data = out_sim,
     type = "l",
     main = "Human prevalence")
plot(z ~ time,
     data = out_sim,
     type = "l",
     main = "Mosquito prevalence")
plot(foi ~ time,
     data = out_sim,
     type = "l",
     main = "Human force of infection")
plot(inf_inc ~ time,
     data = out_sim,
     type = "l",
     main = "Infection incidence")

