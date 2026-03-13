# EpiWave FOI Model

**A Two-Stage Vector-Informed Malaria Transmission Mapping Framework**

## Overview

The **EpiWave FOI Model** is a computational framework for generating malaria risk maps that integrate entomological data (vector abundance, biting rates, mortality) with epidemiological surveillance data. Building on the [epiwave.mapping](https://github.com/idem-lab/epiwave.mapping) framework developed by Nick Golding's group, this project extends it with a mechanistic transmission component (Ross-Macdonald ODE) and vector data integration.

Developed as part of Ernest Moyo's PhD research at NM-AIST in collaboration with Vector Atlas, this framework addresses a critical gap in operational malaria risk mapping: **how to efficiently incorporate vector dynamics into disease models without computational infeasibility**.

### Key Innovation

Traditional mechanistic models face a computational bottleneck when jointly inferring transmission and vector parameters. This project implements a **two-stage framework** that:

1. Uses Vector Atlas data to specify entomological parameters as **fixed inputs**
2. Solves transmission dynamics **once per location** (not at every MCMC iteration)
3. Uses mechanistic predictions as **log-offsets** in a Gaussian Process model
4. Models spatially-correlated residuals between mechanistic prediction and observed data
5. Simultaneously fits to **case counts** (Poisson) and **prevalence surveys** (Binomial)

## Mathematical Framework

### Stage 1: Ross-Macdonald ODE

$$\frac{dx}{dt} = m \cdot a \cdot b \cdot z \cdot (1 - x) - r \cdot x$$

$$\frac{dz}{dt} = a \cdot c \cdot x \cdot (1 - z) - g \cdot z$$

Where $x$ = human infection prevalence, $z$ = mosquito infection prevalence, and parameters $m(t)$, $a(t)$, $g(t)$ are time-varying (from Vector Atlas or temperature models). ITN/IRS intervention effects are applied upstream.

**Mechanistic prediction:**

$$I^*_{s,t} = m_{s,t} \cdot a_{s,t} \cdot b \cdot z_{s,t} \cdot \text{pop}_{s,t}$$

### Stage 2: GP + Dual Likelihood

**Latent infection incidence:**

$$I_{s,t} = \exp\left(\alpha + \log(I^*_{s,t}) + \varepsilon_{s,t}\right)$$

$$\varepsilon \sim \text{GP}(0, K), \quad K = \sigma^2 \cdot K_{\text{space}}(\phi) \cdot K_{\text{time}}(\theta)$$

**Dual likelihood:**

$$C_{s,t} \sim \text{Poisson}(\gamma \cdot I_{s,t}) \quad \text{(case counts)}$$

$$Y_{s,t} \sim \text{Binomial}(N_{s,t},\; x_{s,t}) \quad \text{(prevalence surveys)}$$

**5 free parameters** for MCMC inference:
- $\alpha$ — intercept (log-scale adjustment to mechanistic prediction)
- $\gamma$ — reporting/case ascertainment rate
- $\sigma^2$ — GP marginal variance
- $\phi$ — spatial lengthscale (Matern 5/2)
- $\theta$ — temporal lengthscale (Exponential)

The GP captures spatially-correlated departures from the mechanistic model, while the dual likelihood ensures $\alpha$ and $\gamma$ are individually identifiable.

## Code Architecture

15 functions in `R/epiwave-foi-model.R`:

| Stage | Functions | Purpose |
|-------|-----------|---------|
| **Stage 1** | `get_fixed_m()`, `get_fixed_a()`, `get_fixed_g()` | Entomological parameters |
| | `apply_interventions()` | ITN/IRS effects on m, a, g |
| | `ross_macdonald_ode()` | ODE system definition |
| | `solve_ross_macdonald_multi_site()` | Multi-site ODE solver (deSolve) |
| | `compute_mechanistic_prediction()` | FOI x population -> I* |
| **Stage 2** | `build_gp_kernel()` | Separable Matern 5/2 x Exponential kernel |
| | `simulate_gp_residuals()` | True GP residuals for simulation studies |
| | `simulate_prevalence_surveys()` | Binomial survey data from ODE prevalence |
| | `fit_epiwave_gp()` | GP + dual likelihood via greta/greta.gp |
| **Validation** | `simulate_and_estimate()` | Full pipeline with diagnostic plots |
| | `extract_posterior_summary()` | Posterior processing |
| | `compute_performance_metrics()` | RMSE, coverage, comparison |

## Quick Start

```r
# Install dependencies
install.packages(c("deSolve", "greta", "greta.gp", "MASS", "tidyverse", "ggplot2"))
greta::install_tensorflow()

# Source and run simulation study (start small)
source("R/greta_setup.R")
source("R/epiwave-foi-model.R")
results <- simulate_and_estimate(n_sites = 3, n_times = 12, n_samples = 500)
```

## Project Structure

```
epiwave-foi-model/
├── R/
│   └── epiwave-foi-model.R       # Main implementation
├── papers/
│   └── paper2_epiWaveFOI/        # Objective 2 paper (Rmd + HTML)
├── presentations/
│   └── slidev/                   # Slidev presentation deck
├── correspondence/
│   └── issues/                   # Supervisor feedback threads
├── CLAUDE.md                     # Development constraints
├── PROJECT_DOCUMENTATION.md      # Full technical documentation
├── data/                         # Git-ignored data directory
└── README.md
```

## Documentation

- **[PROJECT_DOCUMENTATION.md](PROJECT_DOCUMENTATION.md)** — Full technical docs (math, API, data requirements)
- **[Paper](papers/paper2_epiWaveFOI/objective2_EpiwaveFOI_Model_Paper.Rmd)** — Objective 2 framework paper

## Author

**Ernest Moyo** — PhD Candidate, NM-AIST / Vector Atlas

- moyoe@nm-aist.ac.tz
- [LinkedIn](https://www.linkedin.com/in/ernest-moyo-96aa3813/)
