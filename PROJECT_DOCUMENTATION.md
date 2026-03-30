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
13. [Model Walkthrough: Baseline Results](#model-walkthrough-baseline-results-2026-03-29)
14. [Appendix](#appendix)

---

# Executive Summary

## Project Overview

The **EpiWave FOI (Force of Infection) Model** is a computational framework for generating high-resolution malaria risk maps that integrate entomological data (vector abundance, biting rates, mortality) with epidemiological surveillance data. Developed as part of PhD research at NM-AIST in collaboration with Vector Atlas, this framework addresses a critical gap in operational malaria risk mapping: **how to efficiently incorporate vector dynamics into geostatistical disease models without computational infeasibility**.

## Key Innovation

Traditional mechanistic disease models face a fundamental computational bottleneck when attempting to jointly infer transmission parameters and vector dynamics from case data. This project implements a **two-stage modeling framework** that:

1. Uses external data sources (Vector Atlas) to specify vector parameters as **fixed inputs** rather than inferring them
2. Solves transmission dynamics **once per location** (not at every MCMC iteration)
3. Uses mechanistic predictions as **log-offsets** in a Gaussian Process model with dual likelihood
4. Models spatially-correlated residuals between mechanistic prediction and observed data via GP
5. Simultaneously fits to **case counts** (Poisson) and **prevalence surveys** (Binomial)

## Operational Impact

- **Spatial mapping:** GP captures where the mechanistic model is wrong, producing corrected risk maps
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
    │ • Vector Atlas    │         │ • GP Residuals       │
    │ • Temperature     │         │ • Dual Likelihood    │
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
| **GP Kernels** | `greta.gp` (R) | Gaussian Process kernels and sparse approximations |
| **Data Manipulation** | `tidyverse` (R) | Data wrangling, tidying, and transformation |
| **Visualization** | `ggplot2` (R) | Publication-quality graphics |
| **GP Simulation** | `MASS` (R) | `mvrnorm()` for simulating true GP residuals |
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
| $Y_{t,s}$ | Prevalence survey positives | Data | $\mathbb{Z}_{\geq 0}$ |
| $N_{t,s}$ | Prevalence survey sample size | Data | $\mathbb{Z}_{> 0}$ |
| $\alpha$ | Intercept (log-scale) | Parameter | $(-\infty, \infty)$ |
| $\gamma$ | Reporting/case ascertainment rate | Parameter | $(0, \infty)$ |
| $\varepsilon_{t,s}$ | GP residual | Latent | $(-\infty, \infty)$ |
| $\sigma^2$ | GP marginal variance | Hyperparameter | $(0, \infty)$ |
| $\phi$ | Spatial lengthscale | Hyperparameter | $(0, \infty)$ |
| $\theta$ | Temporal lengthscale | Hyperparameter | $(0, \infty)$ |

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

## Stage 2: GP + Dual Likelihood

### Model Specification

**Latent infection incidence:**
$$
I_{s,t} = \exp\left(\alpha + \log(I^*_{s,t}) + \varepsilon_{s,t}\right)
$$

Where:
- $\alpha$: Intercept (log-scale systematic adjustment to mechanistic prediction)
- $I^*_{s,t}$: Fixed mechanistic prediction from Stage 1 (used as log-offset)
- $\varepsilon_{s,t}$: Spatially-correlated residual from GP

**GP prior on residuals:**
$$
\varepsilon \sim \text{GP}(0, K), \quad K = \sigma^2 \cdot K_{\text{space}}(\lVert \cdot \rVert; \phi) \cdot K_{\text{time}}(|\cdot|; \theta)
$$

The kernel is **separable**: Matern 5/2 for space × Exponential for time. This matches the epiwave.mapping design and avoids conditioning issues from joint kernels.

### Dual Likelihood

**Case counts (Poisson):**
$$
C_{s,t} \sim \text{Poisson}(\gamma \cdot I_{s,t})
$$

**Prevalence surveys (Binomial):**
$$
Y_{s,t} \sim \text{Binomial}(N_{s,t},\; x_{s,t})
$$

Where $x_{s,t}$ is the ODE human prevalence adjusted by the GP residual.

### Why Dual Likelihood is Required

With case counts only, $\alpha$ (average infection incidence) and $\gamma$ (reporting rate) are **non-identifiable** — they form a ridge in posterior space where $\exp(\alpha) \times \gamma$ is constant. Adding prevalence surveys breaks this because parasite rate surveys are unaffected by case ascertainment, providing independent information about infection incidence.

### Why the GP is Required

The GP models spatially-correlated departures from the mechanistic prediction $I^*$. Without the GP, the model assumes $I^*$ perfectly explains spatial variation — which it does not. The GP captures where and when the mechanistic model is wrong, producing corrected risk maps.

A Negative Binomial likelihood was previously tried (v4) but rejected by the supervisor: it moves spatially-correlated error into the sampling model, preventing spatial mapping of residuals.

### Prior Distributions

| Parameter | Prior | Rationale |
|-----------|-------|-----------|
| $\alpha$ | $\mathcal{N}(0, 1)$ | No systematic bias expected a priori |
| $\gamma$ | $\mathcal{N}(0.1, 0.05)$, truncated > 0 | Literature-informed reporting rate |
| $\log(\sigma^2)$ (`log_sigma`) | $\mathcal{N}(-1, 1)$ | Moderate GP variance |
| $\log(\phi)$ (`log_phi`) | $\mathcal{N}(1, 1)$ | Spatial lengthscale on log scale |
| $\log(\theta)$ (`log_theta`) | $\mathcal{N}(1, 1)$ | Temporal lengthscale on log scale |

### GP Scalability

For computational tractability, sparse GP approximations with inducing points are used via `greta.gp::gp(..., inducing = ...)`. Block-Circulant Embedding (BCB) via FFT is available for larger problems (as in epiwave.mapping).

---

# Implementation Details

## Code Structure

```
epiwave-foi-model/
├── R/
│   ├── epiwave-foi-model.R       # Main implementation (GP + dual likelihood)
│   ├── archive/                   # Previous versions
│   │   └── walkthrough_scripts/  # Step-by-step execution scripts (2026-03-29)
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

#### `fit_epiwave_gp(observed_cases, I_star, x_star, coords, prev_data, ...)`

**Purpose:** Define GP + dual likelihood model with mechanistic offset for Bayesian inference via greta/greta.gp.

> **Note:** The previous function `fit_epiwave_with_offset()` used a NegBin likelihood and has been renamed to `fit_epiwave_negbin()` (kept as deprecated legacy for comparison only).

**Inputs:**
- `observed_cases`: Matrix `[n_times x n_sites]` of case counts
- `I_star`: Matrix `[n_times x n_sites]` of mechanistic predictions
- `x_star`: Matrix `[n_times x n_sites]` of ODE human prevalence
- `coords`: Matrix `[n_obs x 3]` of (lon, lat, time_normalised) coordinates
- `prev_data`: List from `simulate_prevalence_surveys()` (NULL for case-only)
- `use_mechanistic`: Logical — include mechanistic offset?
- `inducing`: Optional matrix of inducing point coordinates for sparse GP

**Output:** `greta` model object (ready for `greta::mcmc()`)

**Model Construction:**
```r
fit_epiwave_gp <- function(observed_cases, I_star, x_star, coords, prev_data, ...) {
  library(greta); library(greta.gp)

  # Priors (5 free parameters)
  alpha     <- normal(0, 1)
  gamma     <- normal(0.1, 0.05, truncation = c(0.001, Inf))
  log_sigma <- normal(-1, 1);  sigma2 <- exp(log_sigma)
  log_phi   <- normal(1, 1);   phi    <- exp(log_phi)
  log_theta <- normal(1, 1);   theta  <- exp(log_theta)

  # Separable GP kernel: Matern 5/2 (space) x Exponential (time)
  K <- build_gp_kernel(sigma2, phi, theta)
  epsilon <- gp(as_data(coords), K, inducing = inducing, tol = 1e-3)

  # Latent incidence with GP residuals
  log_I <- alpha + log(I_star_vec + 1e-10) + epsilon
  I_latent <- exp(log_I)

  # Dual likelihood
  distribution(cases_vec) <- poisson(gamma * I_latent)
  distribution(prev_data$n_positive) <- binomial(prev_data$n_tested, prev_prob)

  model(alpha, gamma, log_sigma, log_phi, log_theta)
}
```

### Simulation-Estimation Study

#### `simulate_and_estimate(n_sites, n_times, run_mcmc, ...)`

**Purpose:** Generate synthetic data with GP residuals and fit competing models for validation.

**Workflow:**
1. **Generate truth:**
   - Known $m, a, g$ with seasonality
   - Apply ITN scale-up
   - Solve ODEs → true $x, z$
   - Compute true $I^*$
   - Simulate true GP residuals $\varepsilon_{\text{true}}$ via `MASS::mvrnorm()`
   - Generate cases $\sim \text{Poisson}(\gamma \cdot \exp(\alpha + \log(I^*) + \varepsilon))$
   - Generate prevalence surveys $\sim \text{Binomial}(N, x \cdot \exp(\varepsilon))$

2. **Fit models (three-way comparison):**
   - **GP+offset:** `fit_epiwave_gp(cases, I_star, ..., use_mechanistic = TRUE)`
   - **GP-only:** `fit_epiwave_gp(cases, I_star, ..., use_mechanistic = FALSE)`
   - **NegBin (legacy):** `fit_epiwave_negbin(cases, I_star)` for comparison

3. **Run MCMC** (if `run_mcmc = TRUE`):
   - 2 chains, 1000 samples each, 500 warmup
   - Extract posterior summaries and performance metrics

4. **Generate 5 diagnostic plots:**
   - p1: Stage 1 vs true incidence (showing gap = GP residuals)
   - p2: True GP residual heatmap (site x time)
   - p3: Alpha-gamma joint posterior (cluster vs ridge)
   - p4: GP hyperparameter posteriors (sigma, phi, theta)
   - p5: Posterior predictive check (GP+offset vs NegBin)

#### `extract_posterior_summary(draws)` and `compute_performance_metrics(draws, true_values)`

**Purpose:** Post-processing of MCMC output for reporting and comparison.

### Validation Results

> **Note:** The results below are from the previous NegBin model (v4) and are retained for historical reference. Updated GP + dual likelihood results will replace these once the new model is validated.

**Previous NegBin results (10 sites, 48 months) — HISTORICAL:**

| Metric | WITH Offset | WITHOUT Offset |
|--------|-------------|----------------|
| **RMSE** | 0.070 | 0.586 |
| **RMSE Improvement** | **88%** | — |
| **Reporting rate** | 0.1003 (true: 0.10) | — |

**GP + dual likelihood validation:** Pending — three-way comparison (GP+offset vs GP-only vs NegBin) in progress.

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

# Step 5: Fit GP model with mechanistic offset and dual likelihood
# model_gp <- fit_epiwave_gp(observed_cases, I_star, x_star = ode_solution$x,
#                             coords = coords, prev_data = prev_data,
#                             use_mechanistic = TRUE, inducing = inducing)
#
# draws_gp <- greta::mcmc(model_gp, n_samples = 1000, warmup = 500, chains = 2)
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

### Stage 2: GP + Dual Likelihood

**Per-iteration complexity:** Depends on GP approximation:
- **Full GP:** $O(N^3)$ — Cholesky decomposition of $N \times N$ covariance matrix
- **Sparse GP (inducing points):** $O(N \cdot M^2)$ where $M$ = number of inducing points
- **BCB (FFT):** $O(N \log N)$ — linear scaling (as in epiwave.mapping)

Parameters: 5 free (alpha, gamma, log_sigma, log_phi, log_theta) + GP latent variables

**Total MCMC time:** Slower than the previous NegBin model due to GP overhead. Depends on:
- Number of observations $N = S \times T$
- GP approximation method and number of inducing points
- Number of iterations and chains

**Benchmark (sparse GP, ~25 inducing points):**
- 3 sites x 12 months, 2 chains x 500 samples: ~1-3 minutes
- 10 sites x 48 months, 2 chains x 1000 samples: ~5-20 minutes

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

**Current:** GP residuals with mechanistic offset

**Extension:** Add environmental covariates $X$ to the GP mean function

**Modified linear predictor:**
$$
\log(I_{s,t}) = \alpha + \log I^*_{s,t} + X_{s,t}^\top \beta + \varepsilon_{s,t}
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

### GP Model Diagnostics

**Check:**
1. α posterior: does the intercept centre near expected log-scale adjustment?
2. γ posterior: does reporting rate concentrate near true value (identifiable via prevalence data)?
3. GP variance σ²: small = mechanistic model explains most variation; large = significant residual structure
4. Lengthscales φ, θ: are spatial/temporal correlation ranges plausible?
5. GP residual surface: does it show interpretable spatial structure?

**Interpretation:**
- WITH offset: GP residuals should be small (σ² near 0) → mechanistic model captures most variation
- WITHOUT offset: GP must absorb all spatial structure → larger σ², wider credible intervals

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
alpha     <- normal(0, 1)           # Log-scale intercept
gamma     <- normal(0.1, 0.05, truncation = c(0.001, Inf))  # Reporting rate
log_sigma <- normal(-1, 1)          # Log GP variance
log_phi   <- normal(1, 1)           # Log spatial lengthscale
log_theta <- normal(1, 1)           # Log temporal lengthscale

# If ESS is low, try more samples
draws <- mcmc(model, n_samples = 2000, warmup = 1000, chains = 4)

# Check: if GP variance posterior is very small, mechanistic model fits well
# This is expected when I* accurately captures spatial variation
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

# Model Walkthrough: Baseline Results (2026-03-29)

This section documents the step-by-step execution of the full pipeline (Stage 1 + Stage 2 MCMC) and the diagnostic findings. It serves as a baseline before any convergence improvements.

## Stage 1 Outputs

Configuration: 10 sites, 49 time points (0–1440 days, monthly), ITN scale-up 0–70%, resistance index 0.2.

| Quantity | Range | Notes |
|----------|-------|-------|
| m (mosquito ratio) | 0.42 – 3.16 | Seasonal + ITN-reduced |
| a (biting rate) | 0.25 – 0.30 | ITN feeding inhibition |
| g (mortality rate) | 0.10 – 0.117 | ITN mortality boost |
| x (human prevalence) | 0.01 – 0.765 | ODE equilibrium |
| z (mosquito prevalence) | 0.001 – 0.640 | ODE equilibrium |
| I* (mechanistic rate) | 0.0004 – 0.468 | m × a × b × z (RATE, not count) |

Stage 1 solves in ~0.8 seconds for 10 sites. ODE is solved once — not at every MCMC iteration.

## Stage 1 Diagnostic Plot: Two-Panel Design

The original single-panel plot mixed rates (I\*, ~0–0.47) with case counts (~0–1850) on the same y-axis, making I\* invisible. Replaced with a two-panel layout:

- **Top panel (Rate Space):** I\* (mechanistic) vs I_true (with GP residuals) on the same scale. Shows exactly what the GP does — exp(epsilon) distorts I\* by 0.14× to 12×.
- **Bottom panel (Count Space):** Expected cases (gamma × I_true × N) vs observed (Poisson draws). Shows the data Stage 2 actually fits to.

This separates the modelling question (can the GP recover the gap?) from the data (what does the model see?).

## GP Residual Structure (True)

True GP parameters: sigma = 0.6, phi = 3.0, rho = 0.75.

- epsilon range: -1.99 to 2.50, mean ≈ 0, SD ≈ 1.02
- Heatmap shows clear spatial correlation (nearby sites share colours) and temporal persistence (AR(1) with rho=0.75)
- 147 prevalence surveys generated (30% of 490 site-time pairs)

## Stage 2 MCMC Results (Baseline)

Settings: 1000 samples, 500 warmup, 2 chains, sparse GP with 25 inducing points.

### GP+offset model (<1% bad transitions, 52 seconds)

| Parameter | True | Posterior Median | Status |
|-----------|------|-----------------|--------|
| alpha | 0.0 | -0.663 | Biased low |
| gamma_rr | 0.1 | 0.153 | Biased high |
| sigma² | 0.36 | 0.789 | Overestimated |
| phi | 3.0 | 1.533 | Underestimated |
| theta (rho) | 0.75 | 0.727 | Recovered well |

### GP-only model (11% bad transitions during warmup, 53 seconds)

| Parameter | Posterior Median | vs GP+offset |
|-----------|-----------------|--------------|
| alpha | -2.159 | Much lower (no offset anchor) |
| gamma_rr | 0.119 | Similar |
| sigma² | 3.685 | ~5× larger — GP absorbs all spatial variation |
| phi | 1.244 | Shorter lengthscale |
| theta | 0.613 | Lower temporal correlation |

### Key Diagnostic Findings

1. **Alpha-gamma ridge (Plot 3):** Clear negative correlation between alpha and gamma in the joint posterior. The dual likelihood (prevalence data) should break this ridge, but 147 surveys across 490 site-time combinations may be insufficient. The posterior does not reach the true values (alpha=0, gamma=0.1).

2. **Multimodal GP hyperparameters (Plot 4):** phi concentrates around 1.2–1.5 but the true value is 3.0. sigma² and theta show multiple modes. Indicates chains have not fully converged with 1000 samples.

3. **Posterior predictive undershoot (Plot 5):** GP+offset prediction captures the general shape but underestimates peaks, consistent with alpha bias pulling predictions down.

4. **GP-only model confirms the framework:** Without the mechanistic offset, GP variance explodes (0.79 → 3.69) because it must explain everything the ODE would have captured. This is the expected behaviour and demonstrates the value of including I\*.

### Issues to Address

- Alpha-gamma non-identifiability ridge persists despite dual likelihood
- GP hyperparameters (phi, sigma²) not recovering true values
- Multimodal posteriors suggest insufficient sampling or poorly-tuned priors
- Need more warmup/samples or structural changes to achieve convergence

Output plots saved in `outputs/plot1_twopanel.png` through `outputs/plot5_posterior_pred.png`.

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
| $\alpha$ | `alpha` | $\mathcal{N}(0, 1)$ | Log-scale intercept adjusting mechanistic prediction |
| $\gamma$ | `gamma` | $\mathcal{N}(0.1, 0.05)$, truncated $> 0$ | Reporting/case ascertainment rate |
| $\log(\sigma^2)$ | `log_sigma` | $\mathcal{N}(-1, 1)$ | GP marginal variance (log scale) |
| $\log(\phi)$ | `log_phi` | $\mathcal{N}(1, 1)$ | Spatial lengthscale, Matérn 5/2 |
| $\log(\theta)$ | `log_theta` | $\mathcal{N}(1, 1)$ | Temporal lengthscale, Exponential |

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

The current GP + dual likelihood model is the result of iterative development:

| Version | Approach | Outcome |
|---------|----------|---------|
| v1 | Joint GP (single RBF kernel) | 490×490 Cholesky near-singular, 0% acceptance |
| v2 | Separable GP (additive kernels) | Ill-conditioned, slow mixing |
| v3 | Hierarchical random effects | α/γ non-identifiable (case-only data) |
| v4 | NegBin with offset (2 params) | Converged but cannot model spatial residuals |
| **v5** | **GP + dual likelihood (5 params)** | **Current — separable kernel, sparse GP, Poisson+Binomial** |

Key lessons incorporated into v5: sparse GP with inducing points (avoids v1/v2 Cholesky issues), dual likelihood with prevalence data (resolves v3 α/γ identifiability), separable Matérn 5/2 × Exponential kernel (following epiwave.mapping). The NegBin v4 was rejected by supervisor because it cannot capture spatially-correlated departures from the mechanistic model.

## D. Glossary

| Term | Definition |
|------|------------|
| **Entomological inoculation rate (EIR)** | Number of infectious mosquito bites per person per unit time |
| **Force of infection (FOI)** | Per-capita rate of infection in susceptible population |
| **Gaussian Process (GP)** | Non-parametric prior over functions; models spatially-correlated residuals |
| **Negative Binomial (NegBin)** | Count distribution allowing overdispersion relative to Poisson (used in deprecated v4) |
| **HMC (Hamiltonian Monte Carlo)** | MCMC algorithm using gradient information |
| **Identifiability** | Whether model parameters can be uniquely determined from data |
| **Offset** | Fixed predictor term in regression model (coefficient = 1) |
| **Ross-Macdonald model** | Compartmental model of vector-borne disease transmission |
| **Vectorial capacity** | Measure of vector population's ability to transmit pathogens |

---

**End of Documentation**

For questions or contributions, contact:
**Ernest Moyo** | moyoe@nm-aist.ac.tz | [GitHub](https://github.com/ernestmoyo/epiwave-foi-model)
