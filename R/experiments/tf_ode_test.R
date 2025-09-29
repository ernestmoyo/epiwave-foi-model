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
# work around a bug stopping the session from working when offline
Sys.setenv(UV_OFFLINE = 1)
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
  plogis(-2 + 0.1 * s)
}

g <- function(day) {
  # model average mosquito lifespan based on seasonality
  s <- seasonality(day)
  lifespan <- exp(2 + 0.1 * s)
  # instantaeous death rate is 1 divided by lifespan
  1 / lifespan
}

# plot these to check

n_years <- 3
day <- seq(0, n_years * 365, by = 30)

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
days_burnin <- 365 * 1
days_sim <- 365 * 4
days_total <- days_burnin + days_sim
days <- seq(0, days_total, by = 30) - days_burnin
init <- c(x = 0.1, z = 0.01)
out <- ode(init, days, deriv, params)

# crop out the burnin
keep <- days > 0
out_sim <- out[keep, ]

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

# 2. Write Tensorflow code for spline interpolation of monthly timeseries data
# to evaluate as a continuous function.

# When fitting the model, we will have at-best monthly estimates of our
# time-varying parameters (mosquito abundance, biting, and survival) or the
# covariates we use to model them. To use computationally-efficient adaptive ODE
# solvers, we need to be able to evaluate those time-varying parameters at any
# time in between these observations, within the tensorflow code and greta
# model. To do this, we can use spline interpolation. We can then write a
# tensorflow function to make computationally efficient evaluations of those
# splines.

# We will use thin-plate splines with radial basis functions, since this option
# allows us to compute the bases for a new time very efficiently inside the
# derivative function. See e.g.:
# https://en.wikipedia.org/wiki/Thin_plate_spline#Radial_basis_function

# Here is example code for thin plate spline interpolation (other types of
# splines can be implemented using e.g. the splines package)

# we observe data on a monthly timestep for two years
months_obs <- 1:24
days_obs <- months_obs * 30
mt_obs <- m(days_obs)

# we want to evaluate on a much smaller timestep
dt <- 0.5
days_pred <- seq(0, 2 * 365, by = dt)

# plot the truth (which we won't have access to in a real application where
# covariates are monthly) and the data we observe as points
par(mfrow = c(1, 1))
plot(m(days_pred) ~ days_pred,
     type = "l")
points(mt_obs ~ days_obs)

# now we build RBF TP splines to estimate the function. First we define the
# knots at which to evaluate the splines, then we compute the basis functions at
# the training times, then we estimate the coefficients using a linear model
# (transformed to have appropriate support), then we compute the basis functions
# at the prediction times, and finally we apply the fitted coefficients to
# estimate the function values at those times

# function to define appropraite knots over a range of values (we could just use
# the months for which we have data, but this is more accurate in general). We
# must pass in the required number of knots, the range of values we want to
# cover (this should be the range of values we want to predict to), and whether
# to place knots on the boundaries (on one each of the lower and upper limits).
define_knots <- function(n_knots, limits, on_boundary = TRUE) {

  # if we are not placing knots on the boundary, add two boundary knots which we
  # will later delete, so that we return n_knots knots
  if(!on_boundary) {
    n_knots <- n_knots + 2
  }

  # space the knots out evenly between the limits
  knots <- seq.int(from = limits[1],
                   to = limits[2],
                   length.out = n_knots)

  # if we are not placing knots on the boundary, remove the two boundary knots
  # so that we return n_knots knots
  if(!on_boundary) {
    knots <- knots[-c(1, n_knots + 2)]
  }
  knots
}

# function to evaluate the basis functions at locations x, given knots, and
# assuming we are working in 1 dimension (the function is different for higher
# dimensions as we must use Euclidean distance)
get_bases <- function(x, knots) {
  diffs <- abs(outer(x, knots, FUN = "-"))
  diffs ^ 2 * log(diffs)
}

# how many knots to use to interpolate this timeseries
n_knots <- length(mt_obs)

# get the knot locations
knots <- define_knots(n_knots, range(days_pred))

# get bases for training data
bases_train <- get_bases(days_obs, knots = knots)

# construct a spline, on the  *log* of m, to ensure that values of m are
# always positive (we will transform back later)
spline_mod <- lm(log(mt_obs) ~ bases_train - 1)

# extract the coefficients
coefs <- spline_mod$coefficients

# get the basis functions at the prediction locations
bases_pred <- get_bases(days_pred, knots = knots)

# evaluate the spline and apply the exponential function to get predictions of
# m, ensuring they are positive
pred <- exp(bases_pred %*% coefs)

lines(pred ~ days_pred,
      col = "red",
      lty = 2,
      lwd = 2)

# This does a good job on this simple and smooth example, even with
# extrapolation at the edges. Though if that chnges with the real data, we can
# fit the spline with additional values on either side (repeating data values if
# needed) to minimise the bias.

# We can simultaneously estimate the coefficients for multiple sites in a
# vectorised manner, and apply them to compute m, by using the linear algebra
# expression for the maximum likelihood solution and the expectation fo the
# linear regression model

# pretend these are the data for N different sites, and that they are different
mt_obs_mat <- cbind(mt_obs, mt_obs, mt_obs)
n_sites <- ncol(mt_obs_mat)

coefs_mat <- solve(t(bases_train) %*% bases_train) %*% t(bases_train) %*% log(mt_obs_mat)

# check the results are very similar
max(abs(coefs_mat[, 1] - coefs))

# predict to values of m for all sites
preds_mat <- exp(bases_pred %*% coefs_mat)

# check this has the right dimensions
identical(dim(preds_mat),
          c(length(days_pred), n_sites))

# plot one of them to confirm it looks right
lines(preds_mat[, 1] ~ days_pred,
      lty = 3,
      lwd = 2,
      col = "blue")


# to use this in practice, we can run most of this code once on the available
# data *before* running the model, to define the knots, training bases, and
# estimate the coefficients for all sites, for each continuous function. E.g.:

#   knots <- define_knots(n_knots, range(days_pred))
#   bases_train <- get_bases(days_obs, knots = knots)
#   coefs_mat <- solve(t(bases_train) %*% bases_train) %*%
#     t(bases_train) %*% log(mt_obs_mat)

# and then *inside the derivative function*, we just need to compute the basis
# functions for the given time point, and evaluate the functions at all sites
# simultaneously, E.g.:

#   bases_t <- get_bases(t, knots = knots)
#   m_t <- exp(bases_t %*% coefs_mat)[1, ]

# e.g.:
t <- 200.67
bases_t <- get_bases(t, knots = knots)
m_t <- exp(bases_t %*% coefs_mat)[1, ]
m_t

# we can write this slightly more efficiently as:
diffs <- abs(knots - t)
bases_t <- diffs ^ 2 * log(diffs)
m_t <- exp(bases_t %*% coefs_mat)[1, ]
m_t

# so we need to convert this bit into Tensorflow code, and write a full
# Tensorflow version of the derivative function

# write a tensorflow function, where t, knots, and coefs all have an initial
# batch dimension (used by greta for vectorisation), and where at each slice of
# the batch dimension, coefs is a matrix with sites as rows and bases as
# columns.

batch_size <- 2L
t_r <- t
knots_r <- knots
coefs_mat_r <- coefs_mat

# make t a matrix, and add a batch dimension
t <- as_tensor(
  array(t_r,
        dim = c(1, 1, 1))
)
# make this a row vector, with a leading dimension for batches
knots <- as_tensor(
  array(knots_r,
        dim = c(1, 1, n_knots))
)
# add initial batch dimension
coefs_mat <- as_tensor(
  array(coefs_mat_r,
        dim = c(1, dim(coefs_mat_r)))
)

# tile most of these by the batch size, as greta will do. This may not be
# needed - it might just tile automatically in the expected way, so long as it
# has the batch dimension is there
knots <- tf$tile(knots, multiples = c(batch_size, 1L, 1L))
coefs_mat <- tf$tile(coefs_mat, multiples = c(batch_size, 1L, 1L))

# get the bases
tf_get_bases <- function(t, knots) {
  diffs <- abs(t - knots)
  diffs ^ 2 * log(diffs)
}

tf_compute_spline <- function(bases, coefs) {
  # compute linear predictor
  pred <- tf$matmul(bases, coefs)
  # transpose back to a column vector
  tf$transpose(pred, perm = c(0L, 2L, 1L))
}

# function to evaluate the spline at a time t, in tensorflow code
tf_evaluate_spline <- function(t, knots, coefs_mat) {
  # compute bases
  bases_t <- tf_get_bases(t, knots)
  tf_compute_spline(bases_t, coefs_mat)
}

s_t <- tf_evaluate_spline(t = t,
                          knots = knots,
                          coefs_mat = coefs_mat)
m_t <- exp(s_t)





# 3. Write Tensorflow code to solve the ODEs with some parameters as continuous
# functions of time.

# function to construct the coefficients for a given observation set and knots,
# where f_observed is an n_sites x n_times matric of observed function values,
# t_observed is an n_times vector of times at which these function values were
# observed, and knots is a vector of times at which spline knots are placed.
get_spline_coefs <- function(f_observed, t_observed, knots) {
  bases_train <- get_bases(t_observed, knots = knots)
  solve(t(bases_train) %*% bases_train) %*% t(bases_train) %*% t(f_observed)
}

# days at which we want to evaluate the integral (monthly over 4y)
days_evaluate <- seq(0, 4*365, by = 30)
# number of days of burnin (1y)
burnin_days <- 365

# number of sites to evaluate
n_sites <- 5

# mock up multisite timeseries data on m*, a*, g* (previously modelled parameter
# values). We should have monthly data within the modelled period, but we could
# replicate data for the burnin period if not
previous_days_observed <- -rev(seq(0, burnin_days, by = 30)[-1])
days_observed <- c(previous_days_observed, days_evaluate)

m_star_observed <- t(replicate(n_sites,
                               m(days_observed)))
a_star_observed <- t(replicate(n_sites,
                               a(days_observed)))
g_star_observed <- t(replicate(n_sites,
                               g(days_observed)))

# define knots across this period to be modelled
n_knots <- length(days_observed)
epsilon <- 1e-3
limits <- c(-burnin_days - epsilon, max(days_evaluate) + epsilon)
knots <- define_knots(n_knots, limits)

# calculate coefficient matrices for these functions, applying appropriate
# transformations to have the right support on the parameters
m_star_coefs <- get_spline_coefs(log(m_star_observed),
                                 days_observed,
                                 knots)
a_star_coefs <- get_spline_coefs(qlogis(a_star_observed),
                                 days_observed,
                                 knots)
g_star_coefs <- get_spline_coefs(log(g_star_observed),
                                 days_observed,
                                 knots)

# greta batch_size



integrate_RMd <- function(
    # burnin period, in the time units of the derivative
    burnin,
    # times to evaluate integral (first value set in the negative () to enable burnin)
    times,
    # initial conditions
    x_0, z_0,
    # tensors (with batch dimensions matching greta) for the model parameters we
    # want to do inference on:
    m_int, m_slope,
    a_int, a_slope,
    g_int, g_slope,
    # scalar r objects for the model parameters we will tret as fixed (we will not
    # evaluate derivatives with respect to these)
    b, c, r,
    # knots for the spline interpolation of continuous functions
    knots,
    # pre-computed coefficient matrices for the predetermined parts of the
    # continuous functions at each site i: m_i(t), a_i(t), g_i(t)
    m_star_coefs, a_star_coefs, g_star_coefs,
    # parameters to control approximation quality (?? - check helpfile when have
    # WiFi)
    rtol = 0.001,
    atol = 1e-06) {

  # transform some things into tensors:

  # convert the fixed parameters to tensors
  b <- as_tensor(b)
  c <- as_tensor(c)
  r <- as_tensor(r)

  # convert the knots and coefficient matrices into tensors, of the correct
  # dimensions: [batch_size, rows, columns]

  # make knots a row vector, with an initial batch dimension, to match matrix
  # equations
  knots <- as_tensor(
    array(knots,
          dim = c(1, 1, length(knots)))
  )

  # add initial batch dimension to each of the coefficient matrices
  m_star_coefs <- as_tensor(
    array(m_star_coefs,
          dim = c(1, dim(m_star_coefs)))
  )

  a_star_coefs <- as_tensor(
    array(a_star_coefs,
          dim = c(1, dim(a_star_coefs)))
  )

  g_star_coefs <- as_tensor(
    array(g_star_coefs,
          dim = c(1, dim(g_star_coefs)))
  )

  # the initial states must have the same dimensions as expected of the
  # derivatives: [batch_size, n_sites, n_states], with n_states = 2 (x then z):
  n_sites <- dim(m_star_coefs)[3]

  # check dims, incl. batch size
  x_0 <- as_tensor(array(rep(x_0, n_sites),
                         dim = c(1, n_sites, 1)))
  z_0 <- as_tensor(array(rep(z_0, n_sites),
                         dim = c(1, n_sites, 1)))

  # combined initial states
  y_0 <- tf$concat(list(x_0, z_0),
                   axis = 2L)

  # convert time parameters
  burnin <- as_tensor(burnin, dtype = tf$int32)
  times <- as_tensor(as_tensor(times, dtype = tf$int32))

  # ODE solver control parameters
  rtol <- as_tensor(rtol)
  atol <- as_tensor(atol)

  # define the derivative function

  # We need to decide how to pass in each parameter object. They only need to be
  # passed in if we want to autodiff them, otherwise we can pass them via lexical
  # scoping. We therefore should not need to pass in the spline objects (knots or
  # each of the coefficient matrices) since these are fixed, not the fixed scalar
  # parameters b, c, and r, but we will need to pass in as arguments the model
  # parameters we want to do inference on: the parameters modifying m(t), a(t),
  # and g(t).

  # For lexical scoping, this function must be defined here, after the
  # parameters that it uses

  tf_deriv <- function(t, y,
                       m_int, m_slope,
                       a_int, a_slope,
                       g_int, g_slope) {

    # t will be dimensionless when used in the ode solver, we need to expand out t
    # to have same dim as a scalar constant (in greta this is 3D) so that it can
    # be used in the same way as the greta array in the R function
    t <- tf$reshape(t, shape = shape(1, 1, 1))

    # evaluate the splines to evaluate the approximated continuous functions of
    # m_star_i(t), a_star_i(t), g_star_i(t), the pre-computed parameters for all
    # sites i at this time
    bases_t <- tf_get_bases(t, knots)
    log_m_star_t <- tf_compute_spline(bases_t, m_star_coefs)
    logit_a_star_t <- tf_compute_spline(bases_t, a_star_coefs)
    log_g_star_t <- tf_compute_spline(bases_t, g_star_coefs)

    # apply the model parameters to adjust these, and then transform each to have
    # the expected support for that parameter (note: tf$nn$sigmoid is the same as
    # the inverse logit function)
    m_t <- exp(m_int + m_slope * log_m_star_t)
    a_t <- tf$nn$sigmoid(a_int + a_slope * logit_a_star_t)
    g_t <- exp(g_int + g_slope * log_g_star_t)

    # pull out x and z from y. These will be in the same orientation as the
    # initial state, so assume this abject has dimensions [batch_size, n_sites,
    # n_states], with n_states = 2 (x then z):
    x <- y[, , 0L, drop = FALSE]
    z <- y[, , 1L, drop = FALSE]

    # derivatives for the Ross MacDonald model of the fractions of
    # people (x) and mosquitoes(z) infected
    dx_dt = m_t * a_t * b * z * (1 - x) - r * x
    dz_dt = a_t * c * x * (1 - z) - g_t * z

    # combine the derivatives along the final dimension to have the same shape as
    # the same as the initial value and the inputs, and return
    tf$concat(list(dx_dt, dz_dt),
              axis = 2L)

  }

  # get tensorflow probability ODE function object
  tf_int <- greta:::tfp$math$ode

  # set up the DoPri solver
  solver <- tf_int$DormandPrince(rtol = rtol, atol = atol)

  # get the solution at these times
  integral <- solver$solve(
    ode_fn = tf_deriv,
    initial_time = -burnin,
    initial_state = y_0,
    solution_times = times,
    constants = list(
      m_int = m_int,
      m_slope = m_slope,
      a_int = a_int,
      a_slope = a_slope,
      g_int = g_int,
      g_slope = g_slope
    ))$states

  # TFP prepends a dimension for the times, so the dimensions of the states are
  # now: [n_times, batch_size, n_sites, n_states].To use this in greta, we must
  # put the batch dimension first. So rearrange to: [batch_size, n_times, n_sites,
  # n_states]
  permutation <- seq_along(dim(integral)) - 1L
  permutation[1:2] <- permutation[2:1]
  integral <- tf$transpose(integral, perm = permutation)
  integral
}

integral <- integrate_RMd(
  # burnin period, in the time units of the derivative
  burnin = burnin_days,
  # times to evaluate integral (first value set in the negative () to enable burnin)
  times = days_evaluate,
  # initial conditions
  x_0 = 0.1,
  z_0 = 0.1,
  # tensors (with batch dimensions matching greta) for the model parameters we
  # want to do inference on:
  m_int = as_tensor(array(0, dim = c(1, 1, 1))),
  m_slope = as_tensor(array(1, dim = c(1, 1, 1))),
  a_int = as_tensor(array(0, dim = c(1, 1, 1))),
  a_slope = as_tensor(array(1, dim = c(1, 1, 1))),
  g_int = as_tensor(array(0, dim = c(1, 1, 1))),
  g_slope = as_tensor(array(1, dim = c(1, 1, 1))),
  # scalar r objects for the model parameters we will tret as fixed (we will not
  # evaluate derivatives with respect to these)
  b = params$b,
  c = params$c,
  r = params$r,
  # knots for the spline interpolation of continuous functions
  knots = knots,
  # pre-computed coefficient matrices for the predetermined parts of the
  # continuous functions at each site i: m_i(t), a_i(t), g_i(t)
  m_star_coefs = m_star_coefs,
  a_star_coefs = a_star_coefs,
  g_star_coefs = g_star_coefs,
  # parameters to control approximation quality (?? - check helpfile when have
  # WiFi)
  rtol = 0.001,
  atol = 1e-06)



# extract timeseries of each state parameter (1st batch, all times, first site,
# then each state)
x_r <- as.array(integral[0L, , 0L, 0L])
z_r <- as.array(integral[0L, , 0L, 1L])

# plot these curves, and overplot the R version with the true (not splined)
# continuous parameters, to visually check the fit

par(mfrow = c(2, 1))
plot(x_r ~ days_evaluate,
     type = "l",
     xlab = "days",
     ylab = "x")
# R version
lines(x ~ time,
     data = out_sim,
     lty = 3,
     lwd = 2,
     col = "red")

plot(z_r ~ days_evaluate,
     type = "l",
     xlab = "days",
     ylab = "x")

# R version
lines(z ~ time,
      data = out_sim,
      lty = 3,
      lwd = 2,
      col = "red")


# do tf_function to compile it
tffun_integrate_RMd <- tensorflow::tf_function(integrate_RMd)
# this needs debugging, and probably an input signature to be efficient!


# make the model parameters have a batch size

# see how the compute time varies with the number of sites

# speed test it, with a batch dimensions and varying number of sites
