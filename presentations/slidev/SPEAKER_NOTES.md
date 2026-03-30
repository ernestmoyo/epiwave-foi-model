# EpiWave FOI Model — Speaker Notes

Framework Development Update, 30 March 2026

---

## SLIDE 1 — Title

Good morning/afternoon everyone. Today I'm presenting a framework development update for my PhD Objective 2 — the EpiWave FOI Model. This is a two-stage vector-informed transmission modelling framework. The goal is to integrate vector spatial data — biting rates, mosquito abundance, survival — into malaria risk mapping approaches, with a case study for Africa, specifically Mozambique or Angola. This is joint work under the supervision of Professor Mirau at NM-AIST, Professor Nick Golding at UWA, Dr David Duncan at Melbourne, Dr Samson Kiware at IHI, and Punam Amratia at the Malaria Atlas Project.

---

## SLIDE 2 — The Challenge

The core problem is that current malaria risk mapping treats vector biology and case data as separate worlds. On one side, we have purely statistical approaches — like MAP's geostatistical models — which fit spatial surfaces to case data but have no mechanistic understanding of how transmission actually works. On the other side, purely mechanistic approaches try to infer all dynamic parameters via MCMC, but that becomes computationally intractable at continental scale. And in between, there's no bridge. Entomological surveillance data — biting rates, mosquito survival, insecticide resistance — collected by vector ecologists remains disconnected from the routine case data that drives risk maps. We need a framework that combines mechanistic understanding with statistical flexibility, without the computational cost.

---

## SLIDE 3 — Our Approach

Our solution is a two-stage framework. Stage 1 is the mechanistic engine. We take fixed entomological parameters from Vector Atlas — mosquito density, biting rate, mortality rate — and solve the Ross-Macdonald ODE system once per spatial unit. This gives us I-star, the mechanistic prediction of expected infection incidence based purely on transmission dynamics. Stage 2 then takes this mechanistic prediction as a fixed offset in a Gaussian Process model. The GP captures spatially-correlated departures between what the mechanistic model predicts and what we actually observe in the data. We use a dual likelihood — Poisson for case counts and Binomial for prevalence surveys — which makes all parameters individually identifiable. The whole Stage 2 has only 5 free parameters for MCMC. This follows the epiwave.mapping framework.

---

## SLIDE 4 — Stage 1: Ross-Macdonald ODE

Here's the mathematical detail of Stage 1. We have the classic Ross-Macdonald coupled ODE system. The first equation describes human infection dynamics — new infections from infectious mosquito bites on susceptible humans, minus recovery. The second equation describes mosquito infection dynamics — mosquitoes becoming infected from biting infected humans, minus mosquito death. The key insight is that we don't infer the entomological parameters. m, a, and g are all fixed from Vector Atlas estimates or temperature-dependent models. We solve these ODEs once per spatial unit to get z — mosquito infection prevalence — then compute the mechanistic prediction I-star equals m times a times b times z. Critically, I-star is a rate, not a count. Population enters the Poisson likelihood later in Stage 2, not here. Interventions — ITNs and IRS — adjust m, a, and g upstream of the ODE through the apply_interventions function, so their effects are automatically reflected in the dynamics.

---

## SLIDE 5 — Stage 1: Code Walkthrough

Let me walk you through the actual R code. Three functions drive Stage 1. The first is ross_macdonald_ode — this is the ODE function that deSolve calls at each integration step. It takes the current state x and z at time t, extracts the time-varying parameters m, a, and g. These can be either constants or interpolation functions — we use R's approxfun to create smooth interpolators from Vector Atlas time series data. Then it computes and returns the derivatives dx/dt and dz/dt.

The second function is solve_ross_macdonald_multi_site. This loops over spatial locations independently. For each site, it creates approxfun interpolators for that site's m, a, and g time series, then calls deSolve's ode function with the LSODA algorithm. Each site is solved once — this is not repeated at every MCMC iteration.

The third function is compute_mechanistic_prediction — a single line. m times a times b times z. That's it. I-star is a rate. Population enters later in the Poisson likelihood. This was a critical design correction — I-star must be a rate, with population entering the Poisson likelihood separately.

---

## SLIDE 6 — Stage 1: Output (Two-Panel Plot)

This is what Stage 1 produces, shown for Site 1 in our simulation study. The top panel shows rate space — comparing I-star, the mechanistic prediction from the ODE, with the true incidence including GP residuals. These are on the same scale, so you can clearly see where and how much the GP distorts the mechanistic prediction. Around months 8 to 12, the GP pushes true incidence well above I-star. Around months 25 to 30, it pushes it below. That gap is the epsilon structure — the spatially-correlated residuals that Stage 2 needs to recover.

The bottom panel shows count space — what the model actually observes. The blue line is expected cases: gamma times I-true times population, which is 0.1 times the rate times 10,000. The black dots are observed cases — Poisson draws around that expectation. You can see the overall declining trend from the simulated ITN scale-up, going from 0% to 70% coverage.

We designed this two-panel view because the original single-panel plot mixed rates and counts on the same axis, which made the mechanistic prediction invisible.

---

## SLIDE 7 — Stage 2: GP + Dual Likelihood

Now Stage 2. The structural equation is log of I equals alpha plus log of I-star plus epsilon. Alpha is an intercept — if the mechanistic model is perfectly calibrated, alpha should be zero. I-star is the fixed offset from Stage 1. And epsilon is the GP residual that captures spatially-correlated departures.

The GP uses a Matern 5/2 spatial kernel with innovations drawn independently at each time step, then correlated across time via an AR(1) process. This matches epiwave.mapping exactly and keeps the latent space at n_sites — about 10 variables — rather than n_sites times n_times which would be 490. That's what makes HMC tractable.

The dual likelihood is critical. Cases follow a Poisson distribution with rate gamma times I times N — gamma is the reporting rate, and population enters here, not in I-star. Prevalence surveys follow a Binomial distribution and are completely independent of gamma. So prevalence data pins down alpha — average infection incidence — while case data informs the product of alpha and gamma. Together, they should make both individually identifiable. I say "should" because we have findings on this that I'll show you shortly.

---

## SLIDE 8 — Stage 2: Code Walkthrough

Here's the actual fit_epiwave_gp function — the heart of Stage 2. At the top, we define five greta variables with their priors. Alpha gets a normal(0,1). Gamma gets a normal centred at 0.1 with standard deviation 0.05, truncated above zero. Sigma-squared and phi are log-normal to ensure positivity. Theta is uniform between 0 and 1 for the AR(1) correlation. All of these are greta variables — TensorFlow computes the gradients automatically for HMC.

Then the GP structure. build_gp_kernel creates a Matern 5/2 kernel with the variance and lengthscale parameters. gp draws spatial innovations — one set per time step. ar1 then correlates these innovations across time. This is ported directly from epiwave.mapping's R code.

The latent incidence is alpha plus log of I-star plus epsilon, exponentiated. Then the dual likelihood: Poisson for cases with expected counts being gamma times I times population, and Binomial for prevalence surveys. The use_mechanistic flag toggles between the WITH offset model and I-star equals zero — the standard geostatistical model — which is exactly the comparison we designed for the simulation-estimation study.

---

## SLIDE 9 — Full Pipeline: Code Architecture

Zooming out, here's the full pipeline as a dependency graph. 15 functions total across three modules. Stage 1 in blue has seven functions — three parameter generators, the intervention function, the ODE function, the multi-site solver, and the mechanistic prediction. Stage 2 in green has the GP kernel builder, the AR(1) function, the prevalence survey simulator, and the main fitting function. The validation layer in purple has the simulation-estimation orchestrator and the metrics functions. The key data flow is: Vector Atlas parameters feed into Stage 1, producing I-star. I-star flows into Stage 2 along with case data and prevalence surveys. MCMC inference on 5 parameters produces posterior predictions.

---

## SLIDE 10 — Stage 2 Evolution

This design didn't come easily — it took six iterations. Versions 1 through 3 tried various GP formulations. Version 1 tried a joint 490-by-490 GP — the Cholesky decomposition was near-singular. Version 2 used additive separable kernels — still ill-conditioned. Version 3 used hierarchical random effects, which converged but alpha and gamma were non-identifiable with case data alone.

Version 4 was a wrong turn — I dropped the GP entirely and used a Negative Binomial likelihood for computational convenience. My supervisor correctly rejected this, pointing out that it assumes the mechanistic model perfectly explains spatial variation, which it doesn't. The GP is needed to model how reality differs from I-star.

Version 5 tried a 3D product kernel GP — mathematically correct but computationally intractable with 100% bad HMC transitions because the latent space was too large at 250 variables.

Version 6 is where we are now. It follows epiwave.mapping directly — spatial GP innovations with AR(1) temporal correlation. Less than 1% bad HMC transitions. We've completed a baseline simulation-estimation run and have diagnosed specific convergence issues, which I'll present next.

---

## SLIDE 11 — Alpha-Gamma Non-Identifiability

This is one of the most important findings. With case data alone, alpha and gamma are fundamentally non-identifiable. For any constant k, you can shift alpha up by k and multiply gamma by e-to-the-minus-k, and get exactly the same case count likelihood. The product exp(alpha) times gamma is identifiable, but the individual parameters are not.

Our design addresses this with a dual likelihood. Prevalence surveys follow a Binomial distribution that doesn't depend on gamma at all. So prevalence data should pin down alpha independently, breaking the ridge.

However — and this is where I need your input — our current simulation with 147 prevalence surveys across 490 site-time combinations still shows a clear ridge. You can see it in this plot — the joint posterior of alpha and gamma shows a clear negative correlation. The dashed lines mark the true values at alpha equals 0 and gamma equals 0.1, and the posterior doesn't quite reach them.

The question is: is 30% survey coverage sufficient to break this ridge? Or do we need more prevalence data? With real data, how much prevalence information would we typically have? This is a legitimate research question — how much prevalence data does the dual likelihood need to work effectively?

---

## SLIDE 12 — GP Residual Structure

This heatmap shows the true GP residuals — the epsilon values we simulated. Each cell is one site at one time point. Blue means true incidence is higher than the mechanistic prediction, red means lower.

Two patterns are visible. First, spatial correlation — at any given time step, nearby sites tend to share similar colours. This comes from the Matern 5/2 kernel with lengthscale phi equals 3. Second, temporal persistence — colours persist across adjacent months, reflecting the AR(1) temporal correlation with theta equals 0.75.

This is the structure that Stage 2 needs to recover. Without this GP, the model would assume I-star perfectly explains spatial variation — which it clearly does not. The GP models where and when reality departs from the mechanistic prediction.

---

## SLIDE 13 — Validation Results: WITH vs WITHOUT

This is the core comparison. On the left, the WITH model uses I-star as a log-offset. On the right, I-star equals zero — the standard geostatistical approach with no mechanistic information.

The key finding: without the mechanistic offset, the GP variance sigma-squared jumps from 0.79 to 3.69 — nearly five times larger. The GP has to absorb all the spatial variation that I-star would have explained. The WITHOUT model also has 11% bad HMC transitions during warmup, compared to less than 1% for the WITH model.

The WITH model recovers the AR(1) temporal correlation well — theta estimated at 0.73 versus the true 0.75. But there are issues. Alpha is estimated at -0.66 instead of 0. Gamma is 0.15 instead of 0.10. These are trading off along the ridge we just discussed. And phi — the spatial lengthscale — is estimated at 1.53 versus the true value of 3.0.

These are not minor discrepancies and I want to be transparent about them. The framework demonstrates the benefit of including vector information — the WITH model is clearly better — but the parameter recovery needs improvement before this is a convincing simulation-estimation result.

---

## SLIDE 14 — MCMC Diagnostics

Here are the actual diagnostic plots. On the left, the GP hyperparameter posteriors. The red dashed lines are the true values. Phi concentrates around 1.2 to 1.5 but the true value is 3.0. Sigma-squared shows multiple modes. Theta is in the right region.

On the right, the posterior predictive check for Site 1. The GP plus offset prediction captures the general seasonal shape but consistently underestimates the peaks, which is consistent with the alpha bias.

We've diagnosed three specific issues. First, phi equals 3.0 on normalised coordinates in a zero-to-one square means minimum correlation of 86% between all site pairs. The data literally cannot distinguish phi equals 3 from phi equals 10 — they produce nearly identical correlation matrices. The prior, which has a median of 1.65, pulls phi toward 1.5.

Second, 1000 samples with 500 warmup may be insufficient. The multimodal posteriors suggest the chains haven't fully explored the posterior.

Third, 147 prevalence surveys may not be enough to break the alpha-gamma ridge fully.

My questions are: should we increase the sampling budget to 4000 or more before diagnosing further? Should we adjust the simulation's true phi to a value the data can actually distinguish? And how much survey coverage is realistic for the study site?

---

## SLIDE 15 — Why This Framework Wins

Despite the convergence issues we're working through, the framework design has four clear advantages. First, computational efficiency — no MCMC over ODE parameters, one-time forward ODE solve per site, only 5 parameters for HMC. Second, biological realism — we're using real entomological data and mechanistic transmission dynamics, not just statistical surfaces. Third, statistical flexibility — the GP captures spatially-correlated departures from the mechanistic prediction, and the dual likelihood makes parameters identifiable. Fourth, operational feasibility — the modular design means you can update vector data independently from case data, generate intervention counterfactuals by re-running Stage 1 only, and the framework scales to national mapping.

---

## SLIDE 16 — Questions and Next Steps

I have four specific questions for you. First, on the alpha-gamma ridge — the dual likelihood hasn't fully broken it at 30% survey coverage. How much prevalence data is needed, and is this a simulation design issue or something more fundamental?

Second, on the simulation's phi value — true phi equals 3.0 on normalised coordinates creates near-total spatial correlation. Should I lower phi or use raw coordinates so the GP has actual spatial variation to estimate?

Third, on sampling — should I invest in 4000 or more samples before making structural changes, or address the priors first?

Fourth, on real data readiness — what data sources should we prioritise for the study site? Publicly available prevalence surveys, routine case data, or both?

My immediate plan is to address the convergence issues, then run a sensitivity study varying survey coverage to quantify how much prevalence data breaks the ridge, scope out available data for the study site, and begin misspecification tests — does the offset help even when I-star is somewhat wrong?

---

## SLIDE 17 — Closing

Thank you for your time. I'm happy to take any questions about the framework, the simulation results, the diagnostics, or the planned next steps.
