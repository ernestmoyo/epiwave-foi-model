# Why the Stage 2 priors look the way they do

This text was moved out of the code comments in `fit_epiwave_gp()` on 2026-10-05. The priors themselves are the `PRIORS` list at the top of `R/epiwave-foi-model.R`.

## The sigma2 / theta ridge (reparameterised 12 Aug 2026, commit `ed7d9a1`)

Convergence under the new priors has not yet been shown at study scale. The R-hat values of 2–4 come from the earlier study.

The GP shape parameters would not converge: R-hat was 2–4. The Workbench check traced this to a variance ridge.

The AR(1) marginal variance is about sigma2 / (1 − theta²). Both the innovation variance sigma2 and the correlation theta drive the one quantity the data can see, so they trade off along a ridge. theta also had a flat prior (`variable(0, 1)`), so nothing held it in place along that ridge.

The fix had two parts:

1. **Reparameterise.** Sample the marginal variance tau2 directly, because that is what the data inform. Then derive sigma2 = tau2 · (1 − theta²). tau2 no longer trades off with theta.
2. **A proper Beta(2, 2) prior on theta.** It replaces the flat prior, keeps theta away from 0 and 1, and gives the sampler something to hold on to.

sigma2 is still the quantity scored against the truth, so the recovery table keeps its meaning.

## Important: which priors produced the paper's numbers

The 50-replicate study in `outputs/sim_estimation_results.RData` ran on 11 Aug 2026. That was before this change. It used:

- sigma2 ~ LogNormal, sampled directly;
- theta ~ Uniform(0, 1).

From now on the harness records `meta$priors` and `meta$git_sha` in every study object, so this cannot be ambiguous again.
