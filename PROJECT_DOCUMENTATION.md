# EpiWave FOI Model: Complete Project Documentation

**Vector-Informed Malaria Transmission Mapping Framework**

**Author:** Ernest Moyo (NM-AIST / Vector Atlas)

---

## Table of Contents

1. [Executive Summary](#executive-summary)
2. [Project Architecture](#project-architecture)
3. [Mathematical Framework](#mathematical-framework)
4. [Implementation Details](#implementation-details)
5. [Data Requirements](#data-requirements)
6. [Usage Guide](#usage-guide)
7. [Performance Characteristics](#performance-characteristics)
8. [Extension Points](#extension-points)
9. [Validation and Diagnostics](#validation-and-diagnostics)
10. [Reproducibility](#reproducibility)
11. [Troubleshooting](#troubleshooting)
12. [References and Resources](#references-and-resources)
13. [Appendix](#appendix)

---

# Executive Summary

## Project Overview

The **EpiWave FOI (Force of Infection) Model** is a computational framework for generating high-resolution malaria risk maps that integrate entomological data (vector abundance, biting rates, mortality) with epidemiological surveillance data. Developed as part of PhD research at NM-AIST in collaboration with Vector Atlas, this framework addresses a critical gap in operational malaria risk mapping: **how to efficiently incorporate vector dynamics into geostatistical disease models without computational infeasibility**.

## Key Innovation

Traditional mechanistic disease models face a fundamental computational bottleneck when attempting to jointly infer transmission parameters and vector dynamics from case data. This project implements a **two-stage modeling framework** that:

1. Uses external data sources (Vector Atlas) to specify vector parameters as **fixed inputs** rather than inferring them
2. Solves transmission dynamics **once per location** (not at every MCMC iteration)
3. Uses mechanistic predictions as **offsets** in a Negative Binomial likelihood model
4. Achieves **17x computational speedup** while maintaining interpretability and flexibility

## Operational Impact

- **Scalability:** National-scale mapping feasible in 1-2 hours (vs. 35+ hours for joint inference)
- **Interpretability:** Parameters retain biological meaning (m = mosquitoes, not abstract coefficients)
- **Flexibility:** Can generate intervention counterfactuals by re-running Stage 1 only
- **Modularity:** Update vector data without re-running statistical inference

---

# Project Architecture

## System Design

```
┌─────────────────────────────────────────────────────────────────┐
│                    EpiWave FOI Model System                      │
└─────────────────────────────────────────────────────────────────┘
                              │
              ┌───────────────┴───────────────┐
              │                               │
              ▼                               ▼
    ┌──────────────────┐          ┌──────────────────┐
    │    STAGE 1       │          │    STAGE 2       │
    │  Mechanistic     │          │  Statistical     │
    │  Prediction      │──────►   │  Inference       │
    └──────────────────┘          └──────────────────┘
              │                               │
              │                               │
    ┌─────────┴─────────┐         ┌──────────┴──────────┐
    │ • Vector Atlas    │         │ • Offset Model      │
    │ • Temperature     │         │ • NegBin Likelihood  │
    │ • Interventions   │         │ • MCMC (HMC)        │
    │ • ODE Solver      │         │ • Uncertainty       │
    └───────────────────┘         └─────────────────────┘
              │                               │
              └───────────────┬───────────────┘
                              ▼
                    ┌─────────────────────┐
                    │  Policy Outputs     │
                    │ • FOI Maps          │
                    │ • Uncertainty       │
                    │ • Counterfactuals   │
                    └─────────────────────┘
```

## Technology Stack

| Component | Technology | Purpose |
|-----------|------------|---------|
| **ODE Solving** | `deSolve` (R) | Fast, adaptive Runge-Kutta solver for Ross-Macdonald equations |
| **Bayesian Inference** | `greta` (R) | Bayesian workflow with TensorFlow backend |
| **Data Manipulation** | `tidyverse` (R) | Data wrangling, tidying, and transformation |
| **Visualization** | `ggplot2` (R) | Publication-quality graphics |
| **Version Control** | Git/GitHub | Code management and collaboration |

## Design Principles

1. **Separation of Concerns:** Mechanistic model (Stage 1) decoupled from statistical model (Stage 2)
2. **Computational Efficiency:** Pre-computation of expensive operations
3. **Reproducibility:** Deterministic workflows with version-controlled code
4. **Extensibility:** Modular design allows easy addition of new features
5. **Operationality:** Designed for use by national malaria control programs

---

# Mathematical Framework

## Notation

| Symbol | Description | Type | Range |
|--------|-------------|------|-------|
| $t$ | Time | Continuous | $[0, T]$ |
| $s$ | Spatial location | Discrete | $s \in \{1, \ldots, S\}$ |
| $x_{t,s}$ | Human infection prevalence | State | $[0, 1]$ |
| $z_{t,s}$ | Mosquito infection prevalence | State | $[0, 1]$ |
| $m_{t,s}$ | Mosquito-to-human ratio | Parameter | $(0, \infty)$ |
| $a_{t,s}$ | Human biting rate (per day) | Parameter | $(0, 1]$ |
| $g_{t,s}$ | Mosquito death rate (per day) | Parameter | $(0, \infty)$ |
| $b$ | Transmission prob (mosquito → human) | Fixed | $[0, 1]$ |
| $c$ | Transmission prob (human → mosquito) | Fixed | $[0, 1]$ |
| $r$ | Human recovery rate (per day) | Fixed | $(0, \infty)$ |
| $I^*_{t,s}$ | Mechanistic incidence prediction | Derived | $[0, \infty)$ |
| $C_{t,s}$ | Observed cases | Data | $\mathbb{Z}_{\geq 0}$ |
| $\lambda$ | Log reporting rate | Parameter | $(-\infty, \infty)$ |
| $\text{size}$ | NegBin overdispersion | Parameter | $(0, \infty)$ |

## Stage 1: Mechanistic Model

### Ross-Macdonald Transmission Dynamics

The foundation is the classic Ross-Macdonald model:

**Human compartment:**
$$
\frac{dx_{t,s}}{dt} = \underbrace{m_{t,s} \cdot a_{t,s} \cdot b \cdot z_{t,s}}_{\text{force of infection}} \cdot (1 - x_{t,s}) - \underbrace{r \cdot x_{t,s}}_{\text{recovery}}
$$

**Mosquito compartment:**
$$
\frac{dz_{t,s}}{dt} = \underbrace{a_{t,s} \cdot c \cdot x_{t,s}}_{\text{acquisition from humans}} \cdot (1 - z_{t,s}) - \underbrace{g_{t,s} \cdot z_{t,s}}_{\text{death}}
$$

**Key Feature:** Time-varying parameters $m_{t,s}$, $a_{t,s}$, $g_{t,s}$ capture:

- **Seasonality** (rainfall-driven mosquito abundance)
- **Interventions** (ITN/IRS effects on vectors)
- **Resistance** (reduced intervention efficacy)

### Force of Infection (FOI)

The FOI represents the per-capita per-day rate at which susceptible humans become infected:

$$
\text{FOI}_{t,s} = m_{t,s} \cdot a_{t,s} \cdot b \cdot z_{t,s}
$$

### Mechanistic Incidence Prediction

Expected number of new infections per time unit:

$$
I^*_{t,s} = \text{FOI}_{t,s} \times \text{pop}_{s,t} = m_{t,s} \cdot a_{t,s} \cdot b \cdot z_{t,s} \cdot \text{pop}_{s,t}
$$

**Implementation Note:** In practice, we work with daily time steps and aggregate to monthly for comparison with surveillance data.

### Intervention Effects

**ITN (Insecticide-Treated Nets):**

Let $n_{t,s} \in [0,1]$ be ITN coverage and $u_{t,s} = 1 - \rho_{t,s} \in [0,1]$ be susceptibility (where $\rho_{t,s}$ is resistance).

$$
\begin{aligned}
m_{t,s}^{\text{ITN}} &= m_{t,s}^{\text{base}} \times \Big[(1-n_{t,s}) + n_{t,s}(1-\pi^m u_{t,s})\Big] \\
a_{t,s}^{\text{ITN}} &= a_{t,s}^{\text{base}} \times \Big[1 - n_{t,s} u_{t,s} \rho^a\Big] \\
g_{t,s}^{\text{ITN}} &= g_{t,s}^{\text{base}} \times \Big[1 + n_{t,s} u_{t,s} \rho^g\Big]
\end{aligned}
$$

Where:
- $\pi^m \in [0,1]$: Proportion of mosquitoes killed/repelled per encounter
- $\rho^a \in [0,1]$: Reduction in successful blood meals
- $\rho^g > 0$: Increase in mortality rate

**Biological Interpretation:**

- **Killing/repelling** ($\pi^m$): ITNs reduce $m$ by removing mosquitoes from population
- **Feeding inhibition** ($\rho^a$): ITNs reduce $a$ by preventing successful bites
- **Contact mortality** ($\rho^g$): ITNs increase $g$ by killing mosquitoes attempting to bite
- **Resistance** ($u$): Reduces all intervention effects proportionally

## Stage 2: Negative Binomial Offset Model

### Model Specification

**Linear predictor:**
$$
\log(\mu_{s,t}) = \lambda + \log(I^*_{s,t})
$$

Where:
- $\lambda$: Log reporting/scaling rate (captures under-reporting, symptomatic fraction)
- $I^*_{s,t}$: Fixed mechanistic prediction from Stage 1 (used as offset)

### Observation Model

**Case counts:**
$$
C_{s,t} \sim \text{NegBin}(\text{size},\; p_{s,t})
$$

Where $p_{s,t} = \frac{\text{size}}{\text{size} + \mu_{s,t}}$.

The Negative Binomial naturally handles overdispersion in case counts. When $\text{size} \to \infty$, the NegBin converges to a Poisson — so the model can recover near-Poisson behaviour when the mechanistic prediction is accurate.

### Parameter Non-Identifiability Resolution

An earlier formulation used separate intercept ($\alpha$) and reporting rate ($\gamma$) parameters:

$$
\log(\mu_{s,t}) = \alpha + \log(\gamma) + \log(I^*_{s,t})
$$

Since $\alpha$ and $\log(\gamma)$ are both additive constants on the log scale, they are **non-identifiable** — only their sum $\lambda = \alpha + \log(\gamma)$ can be estimated from data. The current implementation uses a single merged parameter $\lambda$ (called `log_rate` in code).

### Prior Distributions

| Parameter | Prior | Rationale |
|-----------|-------|-----------|
| $\lambda$ (`log_rate`) | $\mathcal{N}(-2, 1)$ | Centres on ~10% reporting rate; $e^{-2} \approx 0.135$ |
| $\log(\text{size})$ (`log_size`) | $\mathcal{N}(3, 1)$ | Median size ~20; 95% CI from ~1 to ~400 |

**Reparameterisation note:** An earlier version used $\phi \sim \text{Beta}(1, 9)$ for overdispersion, which caused boundary sampling issues (ESS = 82). Reparameterising to $\log(\text{size}) \sim \mathcal{N}(3, 1)$ resolved this (ESS = 489).

---

# Implementation Details

## Code Structure

```
epiwave-foi-model/
├── R/
│   ├── epiwave-foi-model.R       # Main implementation (~580 lines, 11 functions)
│   ├── archive/                   # Previous versions
│   ├── experiments/              # Alternative implementations
│   │   ├── tf_ode_test.R         # TensorFlow ODE solver
│   │   ├── greta_rm_ode_op.R     # Greta ODE operation
│   │   └── speed_test.R          # Performance benchmarking
│   └── ode_speed_test/           # Benchmarking suite
│       ├── ode_speed_test.R
│       └── README.md
├── papers/
│   └── paper2_epiWaveFOI/
│       └── objective2_EpiwaveFOI_Model_Paper.Rmd
├── presentations/
│   └── slidev/                   # Slidev presentation deck
├── data/                         # Git-ignored data directory
├── README.md                     # Project overview
└── .gitignore                    # Excludes data, credentials
```

## Key Functions

### Stage 1: Mechanistic Prediction

#### `get_fixed_m(times, locations, vector_atlas_data = NULL)`

**Purpose:** Obtain mosquito abundance ratios.

**Inputs:**
- `times`: Numeric vector of time points (days)
- `locations`: Character vector of location IDs
- `vector_atlas_data`: Optional data.frame with columns `(time, location, m)`

**Output:** Matrix `[n_times x n_locations]`

**Logic:**
1. If `vector_atlas_data` provided → reshape to matrix
2. Else → generate seasonal pattern: $m(t) = m_0 \times [1 + A \sin(4\pi t/365)]$

**Example:**
```r
times <- seq(0, 365*2, by = 30)  # 2 years, monthly
locations <- c("District_A", "District_B")

# Option 1: Generic seasonality
m_seasonal <- get_fixed_m(times, locations, baseline_m = 2.0, seasonal_amplitude = 0.6)

# Option 2: Vector Atlas data
va_data <- data.frame(
  time = rep(times, length(locations)),
  location = rep(locations, each = length(times)),
  m = rnorm(length(times) * length(locations), mean = 2, sd = 0.5)
)
m_va <- get_fixed_m(times, locations, vector_atlas_data = va_data)
```

#### `apply_interventions(m, a, g, itn_coverage, resistance_index)`

**Purpose:** Modify vector parameters based on intervention coverage.

**Inputs:**
- `m, a, g`: Matrices `[n_times x n_sites]`
- `itn_coverage`: Matrix of ITN coverage `[n_times x n_sites]`, values in `[0,1]`
- `resistance_index`: Scalar or matrix, values in `[0,1]`

**Output:** List `(m = m_adj, a = a_adj, g = g_adj)`

**Example:**
```r
# ITN scale-up: 0% → 70% over 4 years
itn_coverage <- matrix(
  rep(seq(0, 0.7, length.out = n_times), n_sites),
  nrow = n_times, ncol = n_sites
)

# 20% resistance (80% susceptibility)
params_adj <- apply_interventions(
  m = m_baseline,
  a = a_baseline,
  g = g_baseline,
  itn_coverage = itn_coverage,
  resistance_index = 0.2
)
```

#### `solve_ross_macdonald_multi_site(m_matrix, a_matrix, g_matrix, times, ...)`

**Purpose:** Solve ODEs for all sites simultaneously.

**Inputs:**
- `m_matrix, a_matrix, g_matrix`: Parameter matrices `[n_times x n_sites]`
- `times`: Time vector
- `b, c, r`: Transmission parameters (fixed)
- `x0, z0`: Initial conditions

**Output:** List `(x = x_matrix, z = z_matrix)` where each is `[n_times x n_sites]`

**Algorithm:**
```
For each site s:
  1. Create interpolation functions: m(t), a(t), g(t) from matrices
  2. Define ODE system with these time-varying parameters
  3. Solve using deSolve::ode() with LSODA method
  4. Extract x[, s] and z[, s]
```

**Solver Choice:** LSODA is an adaptive solver that switches between stiff and non-stiff methods. Suitable for Ross-Macdonald equations which can become stiff under high transmission or intervention scenarios.

#### `compute_mechanistic_prediction(m_matrix, a_matrix, b, z_matrix, population_matrix)`

**Purpose:** Calculate expected incidence from mechanistic model.

**Formula:** $I^*_{t,s} = m_{t,s} \times a_{t,s} \times b \times z_{t,s} \times \text{pop}_{t,s}$

**Output:** Matrix `[n_times x n_sites]`

### Stage 2: Statistical Inference

#### `fit_epiwave_with_offset(observed_cases, I_star, use_mechanistic = TRUE)`

**Purpose:** Define Negative Binomial model with mechanistic offset for Bayesian inference via greta.

**Inputs:**
- `observed_cases`: Matrix `[n_times x n_sites]` of case counts
- `I_star`: Matrix `[n_times x n_sites]` of mechanistic predictions
- `use_mechanistic`: Logical — include mechanistic offset?

**Output:** `greta` model object (ready for `greta::mcmc()`)

**Model Construction:**
```r
fit_epiwave_with_offset <- function(observed_cases, I_star, use_mechanistic = TRUE) {
  library(greta)
  cases_vec  <- as.vector(observed_cases)
  I_star_vec <- as.vector(I_star)

  # Priors (only 2 free parameters)
  log_rate <- normal(-2, 1)
  log_size <- normal(3, 1)
  size     <- exp(log_size)

  # Linear predictor with mechanistic offset
  if (use_mechanistic) {
    log_mu <- log_rate + log(I_star_vec + 1e-10)
  } else {
    log_mu <- log_rate
  }
  mu   <- exp(log_mu)
  prob <- size / (size + mu)

  # Likelihood
  distribution(cases_vec) <- negative_binomial(size, prob)

  model(log_rate, log_size)
}
```

### Simulation-Estimation Study

#### `simulate_and_estimate(n_sites, n_times, run_mcmc, ...)`

**Purpose:** Generate synthetic data and fit competing models for validation.

**Workflow:**
1. **Generate truth:**
   - Known $m, a, g$ with seasonality
   - Apply ITN scale-up
   - Solve ODEs → true $x, z$
   - Compute true $I^*$
   - Generate cases $\sim \text{NegBin}(\text{size}, p)$ with known reporting rate

2. **Fit models:**
   - **Model A (WITH offset):** `fit_epiwave_with_offset(cases, I_star, use_mechanistic = TRUE)`
   - **Model B (WITHOUT offset):** `fit_epiwave_with_offset(cases, I_star, use_mechanistic = FALSE)`

3. **Run MCMC** (if `run_mcmc = TRUE`):
   - 4 chains, 1000 samples each, 500 warmup
   - Extract posterior summaries and performance metrics

4. **Generate 5 diagnostic plots:**
   - p1: Mechanistic prediction vs observed cases
   - p2: Posterior trace plots
   - p3: WITH vs WITHOUT offset comparison
   - p4: Residual analysis
   - p5: Posterior density plots

#### `extract_posterior_summary(draws)` and `compute_performance_metrics(draws, true_values)`

**Purpose:** Post-processing of MCMC output for reporting and comparison.

### Observed Validation Results

From the simulation-estimation study (10 sites, 48 months):

| Metric | WITH Offset | WITHOUT Offset |
|--------|-------------|----------------|
| **RMSE** | 0.070 | 0.586 |
| **RMSE Improvement** | **88%** | — |
| **Reporting rate** | 0.1003 (true: 0.10) | — |
| **NegBin size** | ~240 (near-Poisson) | ~1.8 (heavy overdispersion) |
| **95% CI width** | 0.0013 | 0.0133 |

**MCMC Convergence (WITH offset model):**

| Parameter | Mean | SD | R-hat | ESS (bulk) |
|-----------|------|-----|-------|------------|
| `log_rate` | -2.30 | 0.020 | 1.004 | 1095 |
| `log_size` | 5.49 | 0.404 | 1.004 | 489 |

---

# Data Requirements

## Stage 1 Inputs

### Vector Data

| Data Type | Source | Format | Temporal Resolution |
|-----------|--------|--------|---------------------|
| Mosquito abundance ($m$) | Vector Atlas, HLC surveys | Raster/shapefile | Monthly preferred |
| Biting rate ($a$) | Literature, temperature models | Scalar or time series | Daily/monthly |
| Death rate ($g$) | Micro-climate models (Mordecai et al.) | Time series | Daily/monthly |

**Vector Atlas Integration:**
- Use VA2 predictions if available
- Fall back to temperature-dependent models for gaps
- Document provenance for reproducibility

### Intervention Data

| Data Type | Source | Format | Notes |
|-----------|--------|--------|-------|
| ITN coverage | DHS, NMCP, MAP | Admin-level time series | Proportion in [0,1] |
| IRS coverage | NMCP, PMI | Admin-level time series | Optional |
| Resistance indices | Bioassay data | Site-specific scalars | From WHO test results |

### Population Data

| Data Type | Source | Format | Notes |
|-----------|--------|--------|-------|
| Population counts | WorldPop, CIESIN | Raster | Annual with projections |
| Age structure | Census, UN | Admin-level | For age-stratified models |

## Stage 2 Inputs

### Epidemiological Data

| Data Type | Source | Format | Quality Requirements |
|-----------|--------|--------|----------------------|
| Monthly cases | DHIS2, NMCP | Admin-level time series | Continuous, minimal missing data |
| Prevalence surveys | DHS, MIS | Point data | Optional, for calibration |

**Data Quality Considerations:**

- **Reporting completeness:** Flag months with known outages
- **Geolocation accuracy:** Exclude facilities with poor GPS precision
- **Diagnostic changes:** Account for shifts (microscopy → RDT)
- **Population denominators:** Ensure consistency with catchment populations

### Spatial Data

| Data Type | Source | Format | Notes |
|-----------|--------|--------|-------|
| Admin boundaries | GADM, HDX | Shapefiles | Match surveillance units |
| Facility locations | MFL, OpenStreetMap | Point data | GPS coordinates |

## Data Preprocessing Pipeline

```r
# 1. Spatial harmonization
admin_grid <- raster::rasterize(admin_polygons, base_raster)

# 2. Temporal alignment
monthly_data <- aggregate_daily_to_monthly(cases_daily)

# 3. Quality control
qc_data <- monthly_data %>%
  filter(reporting_rate > 0.8) %>%  # Exclude low-reporting months
  filter(!is.na(lon), !is.na(lat)) %>%  # Require geolocation
  mutate(cases = pmax(0, cases))  # Non-negative

# 4. Population weighting
pixel_to_admin <- aggregate_population_weighted(
  pixel_data = raster_foi,
  polygons = admin_boundaries,
  population = worldpop_raster
)
```

---

# Usage Guide

## Quick Start

### Installation

```r
# Required packages
install.packages(c("deSolve", "greta", "tidyverse", "ggplot2"))

# TensorFlow backend for greta
greta::install_tensorflow()
```

### Basic Workflow

```r
# Source the model
source("R/epiwave-foi-model.R")

# Step 1: Get vector parameters (example with seasonal pattern)
times <- seq(0, 365*4, by = 30)  # 4 years, monthly
locations <- c("Region_A", "Region_B", "Region_C")

m <- get_fixed_m(times, locations, baseline_m = 2.0, seasonal_amplitude = 0.6)
a <- get_fixed_a(times, locations, baseline_a = 0.3)
g <- get_fixed_g(times, locations, baseline_g = 1/10)

# Step 2: Apply interventions
itn_coverage <- matrix(seq(0, 0.7, length.out = length(times)),
                       nrow = length(times), ncol = length(locations))

params_adj <- apply_interventions(m, a, g, itn_coverage, resistance_index = 0.2)

# Step 3: Solve ODEs
ode_solution <- solve_ross_macdonald_multi_site(
  m_matrix = params_adj$m,
  a_matrix = params_adj$a,
  g_matrix = params_adj$g,
  times = times,
  b = 0.8, c = 0.8, r = 1/7
)

# Step 4: Compute mechanistic prediction
population <- matrix(10000, nrow = length(times), ncol = length(locations))

I_star <- compute_mechanistic_prediction(
  m_matrix = params_adj$m,
  a_matrix = params_adj$a,
  b = 0.8,
  z_matrix = ode_solution$z,
  population_matrix = population
)

# Step 5: Fit NegBin model with mechanistic offset
# model_with <- fit_epiwave_with_offset(observed_cases, I_star, use_mechanistic = TRUE)
# model_without <- fit_epiwave_with_offset(observed_cases, I_star, use_mechanistic = FALSE)
#
# draws_with <- greta::mcmc(model_with, n_samples = 1000, warmup = 500, chains = 4)
# draws_without <- greta::mcmc(model_without, n_samples = 1000, warmup = 500, chains = 4)
```

## Running Simulation Study

```r
# Run full simulation-estimation comparison
results <- simulate_and_estimate(
  n_sites = 10,
  n_times = 48,  # 4 years monthly
  include_interventions = TRUE,
  run_mcmc = TRUE
)

# If run_mcmc = TRUE, results include posterior draws and plots
# If run_mcmc = FALSE, returns Stage 1 output and Plot 1 only

# Extract components
true_incidence <- results$true_incidence
mechanistic_pred <- results$mechanistic_prediction
observed_cases <- results$observed_cases
```

## Generating Counterfactuals

```r
# Scenario 1: Observed (baseline)
I_star_observed <- compute_mechanistic_prediction(...)

# Scenario 2: No interventions
params_no_itn <- list(m = m_baseline, a = a_baseline, g = g_baseline)
ode_no_itn <- solve_ross_macdonald_multi_site(...)
I_star_no_itn <- compute_mechanistic_prediction(...)

# Scenario 3: Universal coverage
itn_universal <- matrix(0.9, nrow = n_times, ncol = n_sites)
params_universal <- apply_interventions(..., itn_coverage = itn_universal)
ode_universal <- solve_ross_macdonald_multi_site(...)
I_star_universal <- compute_mechanistic_prediction(...)

# Compare scenarios
scenario_comparison <- data.frame(
  scenario = rep(c("Observed", "No ITN", "Universal"), each = n_times * n_sites),
  incidence = c(I_star_observed, I_star_no_itn, I_star_universal)
)
```

---

# Performance Characteristics

## Computational Complexity

### Stage 1: Mechanistic Model

**Per-pixel complexity:** $O(T \times N_{\text{ODE steps}})$

- $T$: Simulation duration (days)
- $N_{\text{ODE steps}}$: Adaptive, typically $\sim T$ for LSODA

**Total Stage 1 time:** $O(S \times T)$ where $S$ = number of sites

**Parallelization:** Embarrassingly parallel over sites

**Benchmark:** On Intel i7 (8 cores):
- Single site, 4 years: ~0.5 seconds
- 500 sites, 4 years: ~1 minute (parallelized)

### Stage 2: Negative Binomial Model

**Per-iteration complexity:** $O(N)$ — linear in observations (no covariance matrix)

- $N = S \times T$: Total observations
- Only 2 free parameters (log_rate, log_size) regardless of data size

**Total MCMC time:** Depends on:
- Number of iterations (typically 1,000-4,000)
- Number of chains (typically 4)
- Effective sample size requirements

**Benchmark:**
- 10 sites x 48 months, 4 chains x 1000 samples: ~2-5 minutes
- 100 sites x 48 months, 4 chains x 1000 samples: ~10-15 minutes

## Comparison to Joint Inference

| Approach | ODE Solves | MCMC Parameters | Total Time | Scalability |
|----------|------------|-----------------|------------|-------------|
| **Joint inference** | $N_{\text{iter}} \times S$ | 20+ | 35+ hours | Poor |
| **Two-stage (ours)** | $S$ (once) | **2** | 1-2 hours | Excellent |

---

# Extension Points

## Adding New Features

### 1. Multiple Vector Species

**Current:** Single aggregate vector

**Extension:** Species-specific parameters $m^{(k)}, a^{(k)}, g^{(k)}$ for species $k = 1, \ldots, K$

**Modifications:**
```r
# Stage 1: Solve for each species
z_gambiae <- solve_for_species(m_gambiae, a_gambiae, g_gambiae, ...)
z_funestus <- solve_for_species(m_funestus, a_funestus, g_funestus, ...)

# Combined FOI
FOI_total <- (m_gambiae * a_gambiae * b * z_gambiae) +
             (m_funestus * a_funestus * b * z_funestus)
```

### 2. Age-Stratified Models

**Current:** All-age incidence

**Extension:** Age groups $j = 1, \ldots, J$ with different exposure/susceptibility

**Modifications:**
- $x_{t,s}$ → $x_{j,t,s}$
- Age-specific FOI: $\text{FOI}_{j,t,s} = \psi_j \times \text{FOI}_{t,s}$
- Age-specific case data

### 3. Spatial Covariates in Stage 2

**Current:** NegBin with mechanistic offset only (no spatial residual structure)

**Extension:** Add environmental covariates $X$ and/or spatial random effects

**Modified linear predictor:**
$$
\log(\mu_{s,t}) = \lambda + \log I^*_{s,t} + X_{s,t}^\top \beta
$$

### 4. Non-Stationary Seasonality

**Current:** Fixed seasonal pattern

**Extension:** Year-specific amplitude/phase

**Implementation:** Replace $\sin(2\pi t/365)$ with $\sum_{y} \alpha_y \sin(2\pi t_y / 365)$ where $\alpha_y$ are year-specific

### 5. Temporal Delays

**Current:** Instantaneous reporting

**Extension:** Lag between infection and case presentation

**Implementation:** Convolve $I_t$ with reporting delay distribution $f_\delta$

---

# Validation and Diagnostics

## Model Validation

### 1. Simulation-Estimation Study (Completed)

The primary validation approach generates synthetic data with known parameters and assesses recovery:

- **Stage 1 verification:** Mechanistic prediction I* closely tracks true FOI-driven incidence
- **Parameter recovery:** True reporting rate (0.10) recovered as 0.1003
- **Model comparison:** WITH offset achieves 88% RMSE improvement over WITHOUT offset
- **Uncertainty calibration:** WITH offset produces 10x tighter credible intervals

### 2. Posterior Predictive Checks

**Purpose:** Assess model fit to observed data

```r
# After MCMC, compare posterior predictions to observed
posterior_summary <- extract_posterior_summary(draws)
metrics <- compute_performance_metrics(draws, true_values)
```

### 3. Cross-Validation

**Leave-One-Out (LOO):**
- Hold out one site at a time
- Fit model on remaining sites
- Predict held-out site
- Compute predictive accuracy

**K-Fold Temporal:**
- Split time series into K folds
- Train on past, predict future
- Assess temporal forecasting ability

### 4. Sensitivity Analysis

**Parameters to vary:**
- Fixed transmission parameters ($b, c, r$)
- Prior specifications ($\lambda$, $\text{size}$)
- Vector Atlas uncertainty

**Metrics:**
- Change in posterior means
- Change in interval widths
- Rank correlation of predictions

## Diagnostic Checks

### MCMC Diagnostics

```r
# Convergence: R-hat statistic
rhat <- coda::gelman.diag(draws)
# Should be < 1.01

# Effective sample size
ess <- coda::effectiveSize(draws)
# Should be > 400 per chain

# Trace plots
coda::traceplot(draws)

# Autocorrelation
coda::autocorr.plot(draws)
```

### NegBin Model Diagnostics

**Check:**
1. `log_rate` posterior centres near expected reporting rate?
2. `size` parameter: large = near-Poisson (good mechanistic fit), small = heavy overdispersion
3. Residual patterns across sites and time?

**Interpretation:**
- WITH offset: size ~ 240 (near-Poisson) → mechanistic model captures most variation
- WITHOUT offset: size ~ 1.8 (heavy overdispersion) → model struggles without mechanistic structure

---

# Reproducibility

## Environment Management

```r
# Record session info
session_info <- sessionInfo()
write_rds(session_info, "outputs/session_info.rds")

# Package versions
packrat::snapshot()  # Or renv::snapshot()
```

## Version Control

```bash
# Initialize git
git init
git add R/ papers/ *.Rmd *.R
git commit -m "Initial commit: EpiWave FOI model"

# Remote repository
git remote add origin https://github.com/ernestmoyo/epiwave-foi-model.git
git push -u origin main
```

## Seeds and Reproducibility

```r
# Set seed for all random operations
set.seed(20260202)

# For parallel operations
library(doRNG)
registerDoRNG(seed = 20260202)

# Record in metadata
metadata <- list(
  date = Sys.Date(),
  seed = 20260202,
  r_version = R.version.string,
  package_versions = sapply(loaded_packages, packageVersion)
)
```

---

# Troubleshooting

## Common Issues

### Issue 1: ODE Solver Fails

**Symptoms:** `deSolve::ode()` returns NA or error

**Causes:**
- Parameters out of bounds (negative $m, a, g$)
- Extreme parameter values causing stiffness
- Initial conditions outside [0,1]

**Solutions:**
```r
# Check parameter ranges
stopifnot(all(m_matrix > 0))
stopifnot(all(a_matrix > 0 & a_matrix <= 1))
stopifnot(all(g_matrix > 0))

# Try smaller time step
solution <- ode(..., hmax = 0.1)  # Max step size 0.1 days

# Try different solver
solution <- ode(..., method = "bdf")  # For stiff systems
```

### Issue 2: MCMC Not Converging

**Symptoms:** High R-hat, low ESS, divergent transitions

**Causes:**
- Poor initialization
- Pathological posterior geometry
- Too-vague priors

**Solutions:**
```r
# Current priors (tuned for good convergence)
log_rate <- normal(-2, 1)   # Log reporting rate
log_size <- normal(3, 1)    # Log overdispersion (reparameterised)

# If ESS is low, try more samples
draws <- mcmc(model, n_samples = 2000, warmup = 1000, chains = 4)

# Check: if size posterior is very large (>100), model is near-Poisson
# This is expected when mechanistic prediction is accurate
```

### Issue 3: greta/TensorFlow Setup

**Symptoms:** `greta` fails to initialize or `mcmc()` errors

**Solutions:**
```r
# Reinstall TensorFlow
greta::install_tensorflow()

# Check Python environment
reticulate::py_config()

# Verify greta works
library(greta)
x <- normal(0, 1)
m <- model(x)
draws <- mcmc(m, n_samples = 100)  # Quick test
```

### Issue 4: Stage 2 Model Without Offset Produces Poor Fit

**Expected behaviour:** The WITHOUT-offset model will show:
- Heavy overdispersion (small `size` parameter)
- Wide credible intervals
- Higher RMSE

This is by design — it demonstrates the value of the mechanistic offset. The WITHOUT model serves as a baseline comparison.

---

# References and Resources

## Key Papers

- **Ross-Macdonald model:** Macdonald, G. (1957). The Epidemiology and Control of Malaria. Oxford University Press.
- **Temperature-dependent mosquito traits:** Mordecai et al. (2013). Ecology Letters.
- **ITN effect modeling:** Griffin et al. (2010). PLoS Medicine.
- **Geostatistical disease mapping:** Diggle et al. (1998). JRSS Series C.

## Software Documentation

- **deSolve:** https://cran.r-project.org/web/packages/deSolve/
- **greta:** https://greta-stats.org/
- **Vector Atlas:** https://vectoratlas.icipe.org/

## Data Sources

- **Vector Atlas:** https://vectoratlas.icipe.org/
- **Malaria Atlas Project:** https://malariaatlas.org/
- **WorldPop:** https://www.worldpop.org/
- **DHIS2:** https://dhis2.org/

---

# Appendix

## A. Parameter Values

### Default Fixed Parameters

| Parameter | Symbol | Default Value | Source |
|-----------|--------|---------------|--------|
| Human infectiousness | $b$ | 0.8 | Literature consensus |
| Mosquito competence | $c$ | 0.8 | Laboratory studies |
| Human recovery rate | $r$ | 1/7 per day | ~7 day infectious period |
| Baseline mosquito abundance | $m_0$ | 2.0 | Varies by ecology |
| Baseline biting rate | $a_0$ | 0.3 per day | ~1 bite per 3 days |
| Baseline death rate | $g_0$ | 1/10 per day | ~10 day lifespan |

### ITN Effect Sizes

| Parameter | Symbol | Typical Range | Notes |
|-----------|--------|---------------|-------|
| Kill rate | $\pi^m$ | 0.4-0.7 | From experimental hut trials |
| Feeding inhibition | $\rho^a$ | 0.2-0.4 | Behavioral response |
| Contact mortality | $\rho^g$ | 0.3-0.5 | Additional death risk |

### Stage 2 Priors

| Parameter | Code Name | Prior | Interpretation |
|-----------|-----------|-------|----------------|
| $\lambda$ | `log_rate` | $\mathcal{N}(-2, 1)$ | Log reporting rate; $e^{-2} \approx 0.135$ |
| $\log(\text{size})$ | `log_size` | $\mathcal{N}(3, 1)$ | Overdispersion; median size ~20 |

## B. Computational Requirements

### Minimum Specifications

- **CPU:** Dual-core processor
- **RAM:** 8 GB
- **Storage:** 5 GB for code + intermediate outputs
- **OS:** Windows 10, macOS 10.14+, Ubuntu 18.04+

### Recommended Specifications

- **CPU:** Quad-core or higher
- **RAM:** 16 GB
- **GPU:** Optional, for TensorFlow acceleration
- **Storage:** 50 GB for large-scale analyses

## C. Stage 2 Evolution History

The current NegBin model is the result of iterative development:

| Version | Approach | Outcome |
|---------|----------|---------|
| v1 | Joint GP (5+ params) | Non-identifiable, failed to converge |
| v2 | Separable GP | Better but still non-identifiable |
| v3 | Hierarchical random effects | Acceptance rate ~0% |
| **v4** | **NegBin with offset (2 params)** | **90%+ acceptance, ESS > 400** |

The key insight was that the mechanistic offset already captures the structured variation, so the residual model only needs to handle overdispersion — not spatial-temporal correlation.

## D. Glossary

| Term | Definition |
|------|------------|
| **Entomological inoculation rate (EIR)** | Number of infectious mosquito bites per person per unit time |
| **Force of infection (FOI)** | Per-capita rate of infection in susceptible population |
| **Negative Binomial (NegBin)** | Count distribution allowing overdispersion relative to Poisson |
| **HMC (Hamiltonian Monte Carlo)** | MCMC algorithm using gradient information |
| **Identifiability** | Whether model parameters can be uniquely determined from data |
| **Offset** | Fixed predictor term in regression model (coefficient = 1) |
| **Ross-Macdonald model** | Compartmental model of vector-borne disease transmission |
| **Vectorial capacity** | Measure of vector population's ability to transmit pathogens |

---

**End of Documentation**

For questions or contributions, contact:
**Ernest Moyo** | ernest.moyo@nm-aist.ac.tz | [GitHub](https://github.com/ernestmoyo/epiwave-foi-model)
