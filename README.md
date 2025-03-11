# EpiWave FOI Model

## Overview

The **EpiWave FOI Model** is an epidemiological modeling framework designed to estimate the **Force of Infection (FOI)** in a spatial context using **Bayesian inference** and **mechanistic transmission models**. This script demonstrates:

-   **Synthetic data generation** for spatial locations and malaria transmission covariates.
-   **Ross-Macdonald model** for malaria transmission dynamics.
-   **Mapping spatial covariates** to transmission parameters.
-   **FOI estimation** using mechanistic models and spatial statistical techniques (kriging).
-   **Bayesian inference with Greta** to refine FOI estimates.

This model is an initial **version 1.0** being developed as part of **Ernest Moyo’s PhD research** at the **School of Computational and Communication Science and Engineering** under **NM-AIST** and **Vector Atlas**, focusing on integrating **malaria vector data into risk mapping approaches** in **Tanzania**.

## Authors

-   **Ernest Moyo**
-   **Date:** 2025-03-11

## Dependencies

The script requires the following R packages:

``` r
install.packages(c("deSolve", "dplyr", "sp", "gstat", "ggplot2", "greta", "rootSolve"))
```

Load them in R:

``` r
library(deSolve)    # Solving ODEs
library(dplyr)      # Data wrangling
library(sp)         # Spatial data handling
library(gstat)      # Kriging and spatial statistics
library(ggplot2)    # Data visualization
library(greta)      # Bayesian inference
library(rootSolve)  # Solving equilibrium equations
```

## Workflow

### 1. **Synthetic Data Generation**

-   Generates a grid of **30 spatial locations** in **Tanzania**.
-   Assigns synthetic values for **population, mosquito abundance, and ITN coverage**.
-   Creates synthetic **disease data** (observed cases and prevalence).

### 2. **Ross-Macdonald Transmission Model**

-   Models malaria transmission using a **two-equation ODE system**:

    $$\frac{dx}{dt} = m a b z (1 - x) - r x$$ $$\frac{dz}{dt} = a c x (1 - z) - g z$$

-   Solves the ODE system using **deSolve**.

### 3. **Mapping Spatial Covariates to Transmission Parameters**

-   Assigns transmission parameters ($m, a, r, g$) using **vector index and ITN coverage**.
-   Creates a **spatially heterogeneous transmission model**.

### 4. **Force of Infection (FOI) Estimation**

-   Computes **FOI**: $FOI = m a b z$
-   Uses **rootSolve** to find **equilibrium solutions** for malaria prevalence.
-   Estimates **spatially varying FOI** for each location.

### 5. **Spatial Smoothing with Gaussian Process (Kriging)**

-   Uses **variogram analysis** and **kriging** to smooth FOI estimates.
-   Creates a refined spatial prediction of FOI.

### 6. **Bayesian Inference with Greta**

-   Defines a **hierarchical Bayesian model** linking observed cases to FOI.
-   Uses **MCMC sampling** to infer transmission parameters.

### 7. **Visualization & Outputs**

-   Plots **FOI estimates** across locations.
-   Compares **observed vs. predicted cases**.
-   Generates **spatial maps of FOI** using kriging.

## Future Enhancements

This model is in the **early stages** of development for **Ernest Moyo’s PhD research**, which aims to create a **high-resolution malaria risk mapping framework**. Future enhancements will be closely tied to the objectives of this work:

1.  **Integration with Real-World Vector Atlas Data**
    -   Replace synthetic spatial data with **actual entomological surveillance datasets**.
    -   Incorporate **species-specific vector data** to refine parameter estimates.
2.  **Advanced Spatial Bayesian Modeling**
    -   Develop a **spatiotemporal Bayesian hierarchical model** to capture malaria transmission dynamics more accurately.
    -   Improve uncertainty quantification using **prior knowledge from surveillance programs**.
3.  **Refining Force of Infection Estimation**
    -   Develop a more **mechanistically informed FOI model**, linking **mosquito dynamics, climatic factors, and intervention effects**.
    -   Compare different FOI estimation approaches (e.g., direct ODE solutions vs. Gaussian processes).
4.  **Model Validation & Application in Tanzania**
    -   Use **high-resolution Tanzanian malaria surveillance data** for calibration.
    -   Validate the framework against **field-collected epidemiological data**.
    -   Test intervention **impact scenarios** (e.g., ITN and IRS scale-up).
5.  **Stakeholder Engagement & Policy Integration**
    -   Work with **NMCPs and Vector Control Programs** to integrate findings into operational decision-making.
    -   Develop an **interactive dashboard** to visualize **FOI predictions and risk maps**.
6.  **Scaling Up the Model for Pan-African Use**
    -   Extend the framework beyond Tanzania to **other malaria-endemic African countries**.
    -   Standardize methodologies for **vector-integrated malaria risk mapping**.

## Usage

Run the script in R:

``` r
source("EpiWave_FOI_Model.R")
```

Expected outputs: - Synthetic spatial dataset preview. - ODE solution for the **Ross-Macdonald model**. - Parameter mapping based on covariates. - FOI estimation results. - Bayesian parameter inference summary. - Plots of FOI estimates and spatial distribution.

## Contact

For questions or contributions, please reach out to **Ernest Moyo**.
