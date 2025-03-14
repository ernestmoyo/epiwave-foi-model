# EpiWave FOI Model

## Overview

The **EpiWave FOI Model** provides a structured approach for estimating the **Force of Infection (FOI)** using mechanistic transmission models integrated with Bayesian inference and spatial modeling techniques. This approach facilitates enhanced malaria transmission predictions and informs targeted malaria control strategies.

This document presents the mathematical foundations and workflow clearly without embedding R code, highlighting the model’s theoretical underpinnings aligned with Ernest Moyo’s PhD research at NM-AIST and Vector Atlas, focusing on integrating vector data into malaria risk mapping in Tanzania.

## Mathematical Formulations

### Ross-Macdonald Transmission Dynamics

The Ross-Macdonald model describes malaria transmission between humans and mosquitoes through the following system of ordinary differential equations (ODEs):

$$
\frac{dx}{dt} = m \cdot a \cdot b \cdot z \cdot (1 - x) - r \cdot x
$$

$$
\frac{dz}{dt} = a \cdot c \cdot x \cdot (1 - z) - g \cdot z
$$

Where:

-   $x$ = Proportion of infected humans
-   $z$ = Proportion of infected mosquitoes
-   $m$ = Mosquito abundance relative to humans
-   $a$ = Human biting rate by mosquitoes
-   $b$ = Probability of transmission from mosquito to human per bite
-   $c$ = Probability of transmission from human to mosquito per bite
-   $r$ = Human recovery rate (inverse of infection duration)
-   $g$ = Mosquito death rate (inverse of mosquito lifespan)

### Bayesian Integration of Spatial Covariates

The mosquito abundance parameter $m$ is modeled spatially through Bayesian inference, integrating covariates ($X$):

$$
\log(m) = \alpha + \beta \cdot X
$$

-   $\alpha$ = Intercept parameter (prior: normal distribution)
-   $\beta$ = Covariate effect (positive constraint)

### Force of Infection (FOI) Estimation

The FOI, representing the rate at which susceptible individuals become infected, is estimated as:

$$
\xi = m \cdot a \cdot b \cdot z
$$

A Gaussian Process (GP) is used to model the residual spatial-temporal variability:

$$
\log(\xi) = \log(m \cdot a \cdot b \cdot z) + \epsilon
$$

Where $\epsilon \sim GP(0, K)$, with $K$ being a radial basis function (RBF) kernel.

### Bayesian Hierarchical Model for Observed Cases

Observed malaria cases ($C_{lt}$) in each location are modeled with a Poisson distribution:

$$
C_{lt} \sim \text{Poisson}(\gamma \cdot I_{lt})
$$

Where:

-   $\gamma$ = Symptomatic rate
-   $I_{lt}$ = Total infected individuals estimated through FOI $\xi$

## Workflow Summary

The modeling workflow includes:

1.  **Data Preparation**: Collect and integrate spatial data on vectors, environmental covariates, and malaria incidence.
2.  **Parameter Estimation**: Utilize Bayesian inference to estimate model parameters.
3.  **Dynamic Modeling**: Simulate the Ross-Macdonald ODE system iteratively.
4.  **FOI Calculation**: Estimate the spatial-temporal distribution of FOI incorporating Gaussian Processes.
5.  **Case Simulation and Validation**: Generate predictive distributions for malaria cases and validate using field data.

## Future Enhancements

-   Integration of comprehensive real-world data from: **Vector Atlas, MosquitoDB & NMCP** .
-   Development of advanced spatial-temporal Bayesian hierarchical models.
-   Rigorous validation using Tanzanian epidemiological and entomological datasets.
-   Creation of interactive decision-support dashboards for policy makers.
-   Expansion and standardization of the modeling approach for broader Pan-African implementation.

## Author and Contact

**Ernest Moyo**

-   **Date:** 2025-03-14
-   **Email:** `ernestmoyo35@gmail.com | moyoe@nm-aist.ac.tz`
