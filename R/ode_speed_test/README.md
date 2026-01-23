# ODE Speed Test for EpiWave FOI Model

## Overview

This folder contains scripts for benchmarking different ODE solvers for the Ross-MacDonald malaria transmission model, as requested in [Issue #4](https://github.com/ernestmoyo/epiwave-foi-model/issues/4).

## Objective

Find the most computationally efficient method for solving the Ross-MacDonald ODEs with: - Time-varying parameters (bimodal seasonality) - Multiple sites (N = 10, 50, 100, 500) - 4 years of simulation - Support for automatic differentiation (for Bayesian inference)

## Files

```         
R/ode_speed_test/
├── ode_speed_test.R      # R script using deSolve
├── desolve_results.csv   # Results from deSolve tests
└── speed_comparison.png  # Visualization

python/
└── ode_speed_test.py     # TensorFlow/TFP speed tests
```

## Solvers Tested

### R (deSolve)

-   `lsoda` - Adaptive solver (FORTRAN, switches between stiff/non-stiff)
-   `rk4` - Fixed-step Runge-Kutta 4th order
-   `euler` - Fixed-step Euler method
-   `ode45` - Dormand-Prince adaptive method

### TensorFlow Probability

-   `DormandPrince` - Adaptive Runge-Kutta 4(5)
-   `BDF` - Backward Differentiation Formula (good for stiff systems)
-   `RK4` (custom) - Fixed-step RK4
-   `Euler` (custom) - Fixed-step Euler

## Running the Tests

### R Tests

``` r
source("R/ode_speed_test/ode_speed_test.R")
```

### Python/TensorFlow Tests

``` bash
cd python
pip install tensorflow tensorflow-probability pandas numpy
python ode_speed_test.py
```

## Expected Outcomes

1.  Identify fastest forward-solve method for varying N
2.  Identify fastest method supporting gradient computation
3.  Determine optimal solver for wrapping as greta operation

## Task Checklist (from Issue #4)

-   [x] Set up Ross-MacDonald ODEs with time-varying parameters
-   [x] Write R baseline tests with deSolve\
-   [x] Write TensorFlow code for ODE solving
-   [ ] Speed-test forward-mode solution with different solvers
-   [ ] Speed-test differentiation with different solvers
-   [ ] Wrap optimal solver as greta operation

## Author

Ernest Moyo\
Email: [ernestmoyo35\@gmail.com](mailto:ernestmoyo35@gmail.com){.email}
