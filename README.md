# EpiWave FOI Model

**A Two-Stage Vector-Informed Malaria Transmission Mapping Framework**

## Overview

The **EpiWave FOI Model** is a computational framework for generating malaria risk maps that integrate entomological data (vector abundance, biting rates, mortality) with epidemiological surveillance data. Developed as part of Ernest Moyo's PhD research at NM-AIST in collaboration with Vector Atlas, this framework addresses a critical gap in operational malaria risk mapping: **how to efficiently incorporate vector dynamics into disease models without computational infeasibility**.

### Key Innovation

Traditional mechanistic models face a computational bottleneck when jointly inferring transmission and vector parameters. This project implements a **two-stage framework** that:

1. Uses Vector Atlas data to specify entomological parameters as **fixed inputs**
2. Solves transmission dynamics **once per location** (not at every MCMC iteration)
3. Uses mechanistic predictions as **offsets** in a Negative Binomial likelihood model
4. Achieves **17x speedup** while maintaining biological interpretability

## Mathematical Framework

### Stage 1: Ross-Macdonald ODE

$$\frac{dx}{dt} = m \cdot a \cdot b \cdot z \cdot (1 - x) - r \cdot x$$

$$\frac{dz}{dt} = a \cdot c \cdot x \cdot (1 - z) - g \cdot z$$

Where $x$ = human infection prevalence, $z$ = mosquito infection prevalence, and parameters $m(t)$, $a(t)$, $g(t)$ are time-varying (from Vector Atlas or temperature models). ITN/IRS intervention effects are applied upstream.

**Mechanistic prediction:**

$$I^*_{s,t} = m_{s,t} \cdot a_{s,t} \cdot b \cdot z_{s,t} \cdot \text{pop}_{s,t}$$

### Stage 2: Negative Binomial Calibration

$$\log(\mu_{s,t}) = \lambda + \log(I^*_{s,t})$$

$$C_{s,t} \sim \text{NegBin}(\text{size},\; p_{s,t}), \quad p_{s,t} = \frac{\text{size}}{\text{size} + \mu_{s,t}}$$

Only **2 free parameters** for MCMC inference:
- $\lambda \sim \mathcal{N}(-2, 1)$ — log reporting rate
- $\log(\text{size}) \sim \mathcal{N}(3, 1)$ — overdispersion

## Code Architecture

11 functions in ~580 lines (`R/epiwave-foi-model.R`):

| Stage | Functions | Purpose |
|-------|-----------|---------|
| **Stage 1** | `get_fixed_m()`, `get_fixed_a()`, `get_fixed_g()` | Entomological parameters |
| | `apply_interventions()` | ITN/IRS effects on m, a, g |
| | `ross_macdonald_ode()` | ODE system definition |
| | `solve_ross_macdonald_multi_site()` | Multi-site ODE solver (deSolve) |
| | `compute_mechanistic_prediction()` | FOI × population → I* |
| **Stage 2** | `fit_epiwave_with_offset()` | NegBin model via greta/TensorFlow |
| **Validation** | `simulate_and_estimate()` | Full pipeline with diagnostic plots |
| | `extract_posterior_summary()` | Posterior processing |
| | `compute_performance_metrics()` | RMSE, coverage, comparison |

## Quick Start

```r
# Install dependencies
install.packages(c("deSolve", "greta", "tidyverse", "ggplot2"))
greta::install_tensorflow()

# Source and run simulation study
source("R/epiwave-foi-model.R")
results <- simulate_and_estimate(n_sites = 10, n_times = 48, run_mcmc = TRUE)
```

## Validation Results

From the simulation-estimation study (10 sites, 48 months):

| Metric | WITH Offset | WITHOUT Offset |
|--------|-------------|----------------|
| RMSE | 0.070 | 0.586 |
| **RMSE Improvement** | **88%** | — |
| Reporting rate recovered | 0.1003 (true: 0.10) | — |
| 95% CI width | 0.0013 | 0.0133 |

## Project Structure

```
epiwave-foi-model/
├── R/
│   └── epiwave-foi-model.R       # Main implementation
├── papers/
│   └── paper2_epiWaveFOI/        # Objective 2 paper (Rmd + HTML)
├── presentations/
│   └── slidev/                   # Slidev presentation deck
├── PROJECT_DOCUMENTATION.md      # Full technical documentation
├── data/                         # Git-ignored data directory
└── README.md
```

## Documentation

- **[PROJECT_DOCUMENTATION.md](PROJECT_DOCUMENTATION.md)** — Full technical docs (math, API, data requirements, troubleshooting)
- **[Paper](papers/paper2_epiWaveFOI/objective2_EpiwaveFOI_Model_Paper.Rmd)** — Objective 2 framework paper

## Author

**Ernest Moyo** — PhD Candidate, NM-AIST / Vector Atlas

- ernest.moyo@nm-aist.ac.tz
- [LinkedIn](https://www.linkedin.com/in/ernest-moyo-96aa3813/)
