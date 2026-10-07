// The model code, line by line, beside the maths it computes and a plain note.
// Every code line here must appear verbatim in R/epiwave-foi-model.R; the site
// test (site/tests/walkthrough.test.js) fails if the code and this page drift.
window.EPIWAVE_WALKTHROUGH = [
  {
    id: "constants", fns: ["(top level)"], title: "Fixed constants", file: "Section 1",
    intro: "Values that never change during fitting. b, c and r come from the literature; the time step and the detectability window are modelling choices.",
    rows: [
      { c: "TRANSMISSION <- list(b = 0.5, c = 0.5, r = 1/180, eip_days = 10)", m: "b = c = 0.5,  r = 1/180 day<sup>−1</sup>,  n = 10 days", n: "Transmission probabilities, recovery rate (infections last about six months) and the extrinsic incubation period, all from published studies." },
      { c: "I_STAR_FLOOR <- 1e-6   # keeps log(I*) finite where the ODE gives z ~ 0", m: "I* ≥ 10<sup>−6</sup>", n: "A tiny floor so that log I* never becomes −∞." },
      { c: "DAYS_PER_STEP <- 30    # one model time step (\"month\") = 30 days, for the ODE and the kernel", m: "Δt = 30 days", n: "One model step is a 30-day month." },
      { c: "MAX_DETECT_DAYS <- 30  # how long an infection stays detectable by the test", m: "q(d) = 0 for d > 30", n: "How long after infection a test can still be positive (a placeholder; see the open questions)." }
    ]
  },
  {
    id: "mag", fns: ["get_fixed_m","get_fixed_a","get_fixed_g","m_from_biting_rate"], title: "Mosquito inputs m, a, g", file: "get_fixed_m(), get_fixed_a(), get_fixed_g()",
    intro: "Stand-ins for the Vector Atlas surfaces. Each returns a [sites × months] matrix.",
    rows: [
      { c: "year <- time / 365", m: "t / 365", n: "Convert days to years." },
      { c: "seasonal <- 1 + seasonal_amplitude *\n(sin(4 * pi * year + phase_shift) - 0.3 * cos(2 * pi * year))", m: "s(t) = 1 + A [ sin(4πt/365 + ϕ) − 0.3 cos(2πt/365) ]", n: "Two peaks a year (the sin(4πτ) term) with one stronger season (the cos term): East African rainfall." },
      { c: "m_t <- baseline_m * pmax(seasonal, 0.1)", m: "m(t) = m<sub>0</sub> · max(s(t), 0.1)", n: "Mosquitoes per person through the year, never below 10% of the baseline." },
      { c: "matrix(m_t, nrow = length(location), ncol = length(time), byrow = TRUE)", m: "m<sub>s,t</sub> = m(t) for every site s", n: "The same seasonal curve at every site, one row per site." },
      { c: "matrix(baseline_a, nrow = length(location), ncol = length(time))", m: "a<sub>s,t</sub> = 0.3", n: "Biting rate: bites per mosquito per day (constant stand-in)." },
      { c: "matrix(baseline_g, nrow = length(location), ncol = length(time))", m: "g<sub>s,t</sub> = 0.1", n: "Mosquito death rate per day: a 10-day average life." },
      { c: "biting_rate / a", m: "m = (m·a) / a", n: "The interim Vector Atlas maps give m·a (bites per person per day); dividing by a mapped biting rate a gives m." }
    ]
  },
  {
    id: "interventions", fns: ["itn_effect_retained","apply_interventions"], title: "Bednets and spraying", file: "itn_effect_retained(), apply_interventions()",
    intro: "Interventions change m, a and g before the ODE runs. Resistance lowers the effective net coverage.",
    rows: [
      { c: "1 - 0.46 * (1 - susceptibility)", m: "η(ψ) = 1 − 0.46 (1 − ψ)", n: "Share of a net's effect kept when only a fraction ψ of mosquitoes are susceptible. At ψ = 0 nets still keep 54%, because they block bites." },
      { c: "n <- itn_coverage * itn_effect_retained(itn_susceptibility)", m: "ν = coverage · η(ψ)", n: "Effective net coverage." },
      { c: "m <- m * (1 - n * itn_kill_rate)", m: "m ← m (1 − ν κ)", n: "Nets kill some mosquitoes that try to feed, so fewer mosquitoes per person." },
      { c: "a <- a * (1 - n * itn_feeding_inhibit)", m: "a ← a (1 − ν β)", n: "Nets stop some bites." },
      { c: "g <- g * (1 + n * itn_mortality_boost)", m: "g ← g (1 + ν δ)", n: "Contact with treated nets raises the mosquito death rate." },
      { c: "s <- irs_coverage * irs_susceptibility", m: "ν<sub>IRS</sub> = coverage<sub>IRS</sub> · ψ<sub>IRS</sub>", n: "Effective spraying coverage, with its own susceptibility (IRS insecticides are mostly not pyrethroids)." },
      { c: "g <- g * (1 + s * irs_efficacy)", m: "g ← g (1 + ν<sub>IRS</sub> δ<sub>IRS</sub>)", n: "Sprayed walls kill resting mosquitoes." },
      { c: "a <- a * (1 - s * irs_feeding_inhibit)", m: "a ← a (1 − ν<sub>IRS</sub> β<sub>IRS</sub>)", n: "Spraying slightly reduces biting." }
    ]
  },
  {
    id: "eip", fns: ["eip_gething","eip_survival"], title: "Incubation in the mosquito", file: "eip_gething(), eip_survival()",
    intro: "A mosquito must survive the extrinsic incubation period (EIP) before it can infect anyone.",
    rows: [
      { c: "ifelse(temperature <= 16, Inf, 111 / (temperature - 16))", m: "n(temp) = 111 / (temp − 16),  ∞ for temp ≤ 16 °C", n: "Gething's degree-day rule: the parasite needs 111 degree-days above 16 °C. About 10 days at 27 °C." },
      { c: "(1 + g * eip_days / stages)^(-stages)", m: "S = (1 + g n / k)<sup>−k</sup> ≈ e<sup>−g n</sup>", n: "Share of infected mosquitoes that live through the EIP when it is split into k = 4 stages." }
    ]
  },
  {
    id: "ode", fns: ["ross_macdonald_ode","ross_macdonald_eip_ode"], title: "The Ross-Macdonald equations", file: "ross_macdonald_ode(), ross_macdonald_eip_ode()",
    intro: "x is the share of people infected; z the share of mosquitoes infectious. deSolve calls these functions to get the rates of change.",
    rows: [
      { c: "dx_dt <- m * a * parms$b * z * (1 - x) - parms$r * x", m: "dx/dt = m a b z (1 − x) − r x", n: "People: infectious bites infect the uninfected (1 − x); infected people recover at rate r." },
      { c: "dz_dt <- a * parms$c * x * (1 - z) - g * z", m: "dz/dt = a c x (1 − z) − g z", n: "Mosquitoes (no EIP): biting infected people infects uninfected mosquitoes; mosquitoes die at rate g." },
      { c: "leave <- k * parms$inv_eip(t)   # rate of moving on a stage; 0 when too cold", m: "λ = k / n", n: "With the EIP: each of the k stages is left at rate k/n, so the whole EIP lasts n days on average." },
      { c: "new_infected <- a * parms$c * x * (1 - sum(y) - z)", m: "F = a c x (1 − Σ y<sub>i</sub> − z)", n: "Mosquitoes newly infected per day; only uninfected mosquitoes can be infected." },
      { c: "arriving <- if (i == 1) new_infected else leave * y[i - 1]", m: "in<sub>1</sub> = F,  in<sub>i</sub> = λ y<sub>i−1</sub>", n: "Stage 1 is filled by new infections; each later stage by the stage before it." },
      { c: "dy_dt[i] <- arriving - (g + leave) * y[i]", m: "dy<sub>i</sub>/dt = in<sub>i</sub> − (g + λ) y<sub>i</sub>", n: "A stage loses mosquitoes that die (g) and ones that move on (λ)." },
      { c: "dz_dt <- leave * y[k] - g * z", m: "dz/dt = λ y<sub>k</sub> − g z", n: "Mosquitoes become infectious after the last stage, and keep dying at rate g." }
    ]
  },
  {
    id: "equilibrium", fns: ["rm_equilibrium"], title: "The steady state", file: "rm_equilibrium()",
    intro: "Infections last months, so a simulation must start at the steady state rather than near zero.",
    rows: [
      { c: "S <- eip_survival(g, eip_days, stages)", m: "S = (1 + g n / k)<sup>−k</sup>", n: "Share of infected mosquitoes surviving the EIP (S = 1 when there is no EIP)." },
      { c: "z_of <- function(x) S * a * c * x / (a * c * x + g)", m: "z*(x) = S · a c x / (a c x + g)", n: "With x fixed, the mosquito equations settle at this z." },
      { c: "if (m * a^2 * b * c * S / (g * r) <= 1) {", m: "R₀ = m a² b c S / (g r) ≤ 1", n: "If each infection cannot replace itself, the steady state is no infection." },
      { c: "x <- uniroot(function(x) m * a * b * z_of(x) * (1 - x) - r * x,", m: "solve  m a b z*(x) (1 − x) = r x", n: "Find the human prevalence where new infections balance recoveries." },
      { c: "inflow <- a * c * x * g / (g + a * c * x)", m: "F* = a c x g / (g + a c x)", n: "Daily inflow into the first exposed stage at steady state." },
      { c: "stay <- leave / (leave + g)", m: "λ / (λ + g)", n: "Share of mosquitoes in a stage that move on rather than die." },
      { c: "y <- inflow * stay^(seq_len(stages) - 1) / (leave + g)", m: "y<sub>i</sub>* = F* (λ/(λ+g))<sup>i−1</sup> / (λ+g)", n: "A geometric series: each stage holds a share λ/(λ+g) of the one before." }
    ]
  },
  {
    id: "solver", fns: ["solve_ross_macdonald_multi_site","compute_mechanistic_prediction"], title: "Solving Stage 1 once per site", file: "solve_ross_macdonald_multi_site(), compute_mechanistic_prediction()",
    intro: "No inference happens here: the ODE is solved once, before any sampling.",
    rows: [
      { c: "parms <- list(m = approxfun(times, m_matrix[site, ], rule = 2),", m: "m(t), a(t), g(t) by linear interpolation", n: "Monthly values become smooth functions of time for the ODE solver." },
      { c: "start_state <- rm_equilibrium(m_matrix[site, 1], a_matrix[site, 1], g_matrix[site, 1],", m: "(x<sub>0</sub>, y<sub>0</sub>, z<sub>0</sub>) = steady state at month 0", n: "Start each site at the steady state of its month-0 conditions." },
      { c: "solution <- ode(y = c(x = start_state$x, y = start_state$y, z = start_state$z),", m: "integrate dx/dt, dy/dt, dz/dt from t = 0", n: "deSolve (lsoda) solves the equations over the months." },
      { c: "m_matrix * a_matrix * b * z_matrix", m: "I*<sub>s,t</sub> = m a b z", n: "The force of infection: new infections per susceptible person per day. This is the mechanistic prediction passed to Stage 2." }
    ]
  },
  {
    id: "ar1", fns: ["ar1"], title: "AR(1) through time", file: "ar1() — epiwave.mapping",
    intro: "Links each month's residual to the last. Written as one matrix multiply so TensorFlow runs it fast.",
    rows: [
      { c: "t_mat <- outer(t_seq, t_seq, FUN = \"-\")", m: "D<sub>t,u</sub> = t − u", n: "How many months separate every pair of times." },
      { c: "mask <- lower.tri(t_mat, diag = TRUE)", m: "keep u ≤ t", n: "A month can only be influenced by itself and earlier months." },
      { c: "rho_mat <- (rho ^ t_mat) * mask", m: "R<sub>t,u</sub> = θ<sup>t−u</sup> for u ≤ t", n: "Influence fades by a factor θ for every month back." },
      { c: "t(rho_mat %*% t(innovations))", m: "ε<sub>t</sub> = Σ<sub>u≤t</sub> θ<sup>t−u</sup> f<sub>u</sub>", n: "The same as ε<sub>t</sub> = θ ε<sub>t−1</sub> + f<sub>t</sub>, computed in one step." }
    ]
  },
  {
    id: "detect", fns: ["default_q_daily","build_convolution_matrix"], title: "From infections to test positivity", file: "default_q_daily(), build_convolution_matrix()",
    intro: "Survey prevalence is built from recent infections, weighted by how likely each is still to test positive.",
    rows: [
      { c: "up   <- plogis(days * 2 - 5)", m: "u(d) = logit<sup>−1</sup>(2d − 5)", n: "Detectability rises over the first few days after infection." },
      { c: "down <- 1 - plogis(days / 2 - 10)", m: "v(d) = 1 − logit<sup>−1</sup>(d/2 − 10)", n: "…and fades later." },
      { c: "q <- up * down", m: "q(d) = u(d) · v(d)", n: "Probability a person infected d days ago still tests positive." },
      { c: "if (i >= 1) C[i, j] <- weights[lag + 1]", m: "C<sub>t−ℓ, t</sub> = w<sub>ℓ</sub>", n: "Place each monthly weight w<sub>ℓ</sub> so that (I C)<sub>t</sub> = Σ<sub>ℓ</sub> I<sub>t−ℓ</sub> w<sub>ℓ</sub>. The weights come from Nick's transform_convolution_kernel(), which integrates q(d) over each month." }
    ]
  },
  {
    id: "fit", fns: ["fit_epiwave_gp"], title: "Stage 2: the greta model", file: "fit_epiwave_gp()",
    intro: "Every line from priors to likelihood. This is the statistical model that is fitted by MCMC.",
    rows: [
      { c: "n_times <- ncol(observed_cases)", m: "", n: "Number of months." },
      { c: "alpha    <- normal(PRIORS$alpha$mean, PRIORS$alpha$sd)", m: "α ~ Normal(0, 1)", n: "Overall level of incidence relative to I*." },
      { c: "gamma_rr <- normal(PRIORS$gamma_rr$mean, PRIORS$gamma_rr$sd,", m: "γ ~ Normal(0.1, 0.05), γ > 0.001", n: "Reporting rate: share of infections that appear as reported cases." },
      { c: "tau2     <- lognormal(PRIORS$tau2$meanlog, PRIORS$tau2$sdlog)", m: "τ² ~ LogNormal(−0.5, 0.5)", n: "Overall variance of the residual field, the quantity the data see." },
      { c: "theta    <- beta(PRIORS$theta$shape1, PRIORS$theta$shape2)", m: "θ ~ Beta(2, 2)", n: "Month-to-month persistence of the residuals." },
      { c: "sigma2   <- tau2 * (1 - theta ^ 2)", m: "σ² = τ² (1 − θ²)", n: "The innovation variance follows from τ² and θ; sampling τ² avoids the σ²–θ ridge." },
      { c: "phi      <- lognormal(PRIORS$phi$meanlog, PRIORS$phi$sdlog)", m: "φ ~ LogNormal(0.5, 0.5)", n: "Spatial lengthscale: how far residual similarity reaches." },
      { c: "kernel  <- mat52(lengthscales = phi, variance = sigma2)", m: "K(h) = σ² (1 + √5 h/φ + 5h²/(3φ²)) e<sup>−√5 h/φ</sup>", n: "Matérn 5/2 covariance between two sites a distance h apart." },
      { c: "f       <- gp(spatial_coords, kernel, inducing = inducing, n = n_times, tol = gp_tol)", m: "f<sub>t</sub> ~ GP(0, K) for each month t", n: "A fresh spatial pattern for every month." },
      { c: "epsilon <- ar1(rho = theta, innovations = f)", m: "ε<sub>t</sub> = θ ε<sub>t−1</sub> + f<sub>t</sub>", n: "Chain the monthly patterns through time." },
      { c: "if (center_alpha) epsilon <- epsilon - mean(epsilon)", m: "ε ← ε − mean(ε)", n: "Optional: centre the field so α is the only intercept (α then means α + mean ε)." },
      { c: "log_I <- alpha + log(pmax(I_star, I_STAR_FLOOR)) + epsilon", m: "log I = α + log I* + ε", n: "The core equation: incidence is the mechanistic prediction, scaled by e<sup>α</sup> and adjusted where the mechanism is wrong." },
      { c: "log_I <- alpha + epsilon   # I* = 0: the standard geostatistical model", m: "log I = α + ε", n: "The comparison model with I* = 0: the standard geostatistical model." },
      { c: "I_latent <- exp(log_I)", m: "I = e<sup>log I</sup>", n: "Back to a rate per person per day." },
      { c: "expected_cases <- gamma_rr * I_latent * N_pop * DAYS_PER_STEP", m: "μ = γ · I · N · 30", n: "Expected reported cases per site-month: rate × population × days, times the reporting rate." },
      { c: "cases <- as_data(observed_cases)", m: "C (data)", n: "Hand the observed case counts to greta as data." },
      { c: "distribution(cases) <- poisson(expected_cases)", m: "C ~ Poisson(μ)", n: "Case-count likelihood." },
      { c: "cases <- as_data(observed_cases[observed_sites, , drop = FALSE])\ndistribution(cases) <- poisson(expected_cases[observed_sites, ])", m: "C<sub>s</sub> ~ Poisson(μ<sub>s</sub>), observed sites only", n: "With held-out sites, only observed sites enter the likelihood; the GP still predicts the others." },
      { c: "survey_site <- (prev_data$survey_indices - 1) %% nrow(observed_cases) + 1", m: "site of survey j", n: "Which site each survey belongs to (used to drop held-out sites)." },
      { c: "survey_used <- seq_along(survey_site)\nsurvey_used <- which(survey_site %in% observed_sites)", m: "", n: "Use every survey, or only those at observed sites." },
      { c: "survey_cells <- prev_data$survey_indices[survey_used]", m: "(s, t) of each survey", n: "The site-month of each survey used." },
      { c: "detectable <- I_latent %*% as_data(prev_conv_matrix)", m: "Λ<sub>t</sub> = Σ<sub>ℓ</sub> I<sub>t−ℓ</sub> w<sub>ℓ</sub>", n: "Expected detectable infections per person: recent incidence weighted by detectability." },
      { c: "p_positive <- 1 - exp(-detectable[survey_cells])", m: "p = 1 − e<sup>−Λ</sup>", n: "Chance a tested person is positive: at least one detectable infection. Depends on I, so on α and ε." },
      { c: "positives  <- as_data(prev_data$n_positive[survey_used])", m: "Y (data)", n: "Number testing positive in each survey." },
      { c: "distribution(positives) <- binomial(prev_data$n_tested[survey_used], p_positive)", m: "Y ~ Binomial(T, p)", n: "Survey likelihood. With the cases, it separates α from γ." },
      { c: "greta_model <- model(alpha, gamma_rr, sigma2, phi, theta)", m: "posterior ∝ likelihood × priors", n: "Collect the parameters to sample." },
      { c: "attr(greta_model, \"I_latent\") <- I_latent", m: "", n: "Keep the latent incidence so a fit can be turned into a map." }
    ]
  },
  {
    id: "simulate", fns: ["simulate_epiwave_data","simulate_gp_residuals"], title: "Simulating the truth", file: "simulate_epiwave_data(), simulate_gp_residuals()",
    intro: "How the test data are made, so the fit can be checked against known values.",
    rows: [
      { c: "itn_coverage <- outer(itn_max, seq(0, 1, length.out = length(times)))", m: "ν<sub>s,t</sub> = ν<sub>max,s</sub> · t / n<sub>times</sub>,  t = 0, …, n<sub>times</sub>", n: "Net coverage rises in a straight line to each site's maximum." },
      { c: "I_star <- compute_mechanistic_prediction(m, a, tp$b, ode$z)", m: "I* = m a b z", n: "Stage 1 output for the simulated sites." },
      { c: "K_space <- (1 + r + r^2 / 3) * exp(-r) + diag(1e-6, n_sites)", m: "K<sub>ij</sub> = (1 + r + r²/3) e<sup>−r</sup>,  r = √5 h<sub>ij</sub>/φ", n: "Spatial correlation between the simulated sites: the same Matérn 5/2 as the model, with σ² = 1 (σ is applied later). Here r is the scaled distance, not the recovery rate. A tiny jitter keeps it invertible." },
      { c: "innovations[, t] <- MASS::mvrnorm(1, mu = rep(0, n_sites), Sigma = K_space)", m: "u<sub>t</sub> ~ MVN(0, K)", n: "A correlated random pattern for each month." },
      { c: "epsilon[, t] <- rho * epsilon[, t - 1] + sigma * innovations[, t]", m: "ε<sub>t</sub> = θ ε<sub>t−1</sub> + σ u<sub>t</sub>", n: "The true residuals, built with plain R so the truth does not share code with the model being tested." },
      { c: "I_true <- exp(tp$alpha + log(pmax(I_star, I_STAR_FLOOR)) + epsilon)", m: "I<sup>true</sup> = e<sup>α</sup> I* e<sup>ε</sup>", n: "True incidence." },
      { c: "expected_cases <- tp$reporting_rate * I_true * population * DAYS_PER_STEP", m: "μ = γ · I · N · 30", n: "Expected reported cases." },
      { c: "observed_cases <- matrix(rpois(length(expected_cases), expected_cases),", m: "C ~ Poisson(μ)", n: "Draw the observed case counts." },
      { c: "prevalence_true <- 1 - exp(-(I_true %*% conv_matrix))", m: "p = 1 − e<sup>−Λ</sup>", n: "True survey prevalence, by the same rule the model uses." }
    ]
  }
];
