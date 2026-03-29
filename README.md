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

**Mechanistic prediction (infection incidence RATE):**

$$I^*_{s,t} = m_{s,t} \cdot a_{s,t} \cdot b \cdot z_{s,t}$$

I* is a rate, not a count. Population enters the Poisson likelihood separately.

### Stage 2: Spatial GP + AR(1) Temporal + Dual Likelihood

Following [epiwave.mapping](https://github.com/idem-lab/epiwave.mapping) directly:

**Latent infection incidence:**

$$I_{s,t} = \exp\left(\alpha + \log(I^*_{s,t}) + \varepsilon_{s,t}\right)$$

**Spatial GP innovations + AR(1) temporal correlation:**

$$f_t \sim \text{GP}(0, \sigma^2 \cdot \text{Matérn}_{5/2}(\phi))$$

$$\varepsilon_t = \theta \cdot \varepsilon_{t-1} + f_t$$

**Dual likelihood:**

$$C_{s,t} \sim \text{Poisson}(\gamma \cdot I_{s,t} \cdot N_s) \quad \text{(case counts)}$$

$$Y_{s,t} \sim \text{Binomial}(N_{s,t},\; x_{s,t}) \quad \text{(prevalence surveys)}$$

**5 free parameters** for MCMC inference:
- $\alpha$ — intercept (log-scale adjustment to mechanistic prediction)
- $\gamma$ — reporting/case ascertainment rate
- $\sigma^2$ — GP marginal variance
- $\phi$ — spatial lengthscale (Matern 5/2)
- $\theta$ — AR(1) temporal correlation [0, 1]

The GP captures spatially-correlated departures from the mechanistic model. The dual likelihood ensures $\alpha$ and $\gamma$ are individually identifiable. The simulation-estimation study compares WITH I* offset vs I*=0 (standard geostatistical, no mechanistic information), as specified by Nick Golding.

## Code Architecture

15 functions in `R/epiwave-foi-model.R`:

| Stage | Functions | Purpose |
|-------|-----------|---------|
| **Stage 1** | `get_fixed_m()`, `get_fixed_a()`, `get_fixed_g()` | Entomological parameters |
| | `apply_interventions()` | ITN/IRS effects on m, a, g |
| | `ross_macdonald_ode()` | ODE system definition |
| | `solve_ross_macdonald_multi_site()` | Multi-site ODE solver (deSolve) |
| | `compute_mechanistic_prediction()` | I* = m·a·b·z (incidence rate) |
| **Stage 2** | `build_gp_kernel()` | Spatial Matern 5/2 kernel |
| | `ar1()` | AR(1) temporal correlation (from epiwave.mapping) |
| | `simulate_gp_residuals()` | Spatial GP + AR(1) residuals for simulation |
| | `simulate_prevalence_surveys()` | Binomial survey data from ODE prevalence |
| | `fit_epiwave_gp()` | GP + dual likelihood via greta/greta.gp |
| **Validation** | `simulate_and_estimate()` | WITH offset vs I*=0 comparison |
| | `extract_posterior_summary()` | Posterior processing |
| | `compute_performance_metrics()` | RMSE, coverage, comparison |

## Quick Start

```r
# Install R dependencies
install.packages(c("deSolve", "greta", "greta.gp", "MASS", "dplyr", "tidyr", "ggplot2"))

# Python/TF setup (see greta_setup.R for details)
Sys.setenv(RETICULATE_PYTHON = "path/to/r-greta/python.exe")

# Source and run simulation study
source("R/greta_setup.R")
source("R/epiwave-foi-model.R")
results <- simulate_and_estimate(n_sites = 10, n_times = 24)
```

## Project Structure

```
epiwave-foi-model/
├── R/
│   ├── epiwave-foi-model.R       # Main implementation (Stage 1 + Stage 2)
│   ├── greta_setup.R             # Python/TF/greta environment setup
│   └── experiments/
│       ├── sim_estimation_study.R # Formal simulation-estimation study
│       └── speed_test.R          # ODE solver speed benchmarks
├── docs/
│   └── development_journey.md    # Full v1-v6 iteration history
├── papers/
│   └── paper2_epiWaveFOI/        # Objective 2 paper (Rmd + HTML)
├── presentations/
│   └── slidev/                   # Slidev deck with animations
├── outputs/                      # Saved results, plots (git-ignored)
├── correspondence/
│   └── issues/                   # Supervisor feedback threads
├── meetings/
│   └── PhD_Work/                 # Nick's brainwave, panel notes
└── README.md
```

## Documentation

- **[docs/development_journey.md](docs/development_journey.md)** — Full v1-v6 iteration history with Mermaid diagrams
- **[PROJECT_DOCUMENTATION.md](PROJECT_DOCUMENTATION.md)** — Technical documentation
- **[Paper](papers/paper2_epiWaveFOI/objective2_EpiwaveFOI_Model_Paper.Rmd)** — Objective 2 framework paper

## Author

**Ernest Moyo** — PhD Candidate, NM-AIST / Vector Atlas

- moyoe@nm-aist.ac.tz
- [LinkedIn](https://www.linkedin.com/in/ernest-moyo-96aa3813/)
