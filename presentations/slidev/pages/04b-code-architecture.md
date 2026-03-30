---
layout: default
---

# Code Architecture: `epiwave-foi-model.R`

15 functions — modular, two-stage pipeline with GP residuals and dual likelihood

```mermaid {scale: 0.45}
graph LR
    subgraph S1["STAGE 1: Mechanistic Engine"]
        M["get_fixed_m()"] --> INT["apply_interventions()"]
        A["get_fixed_a()"] --> INT
        G["get_fixed_g()"] --> INT
        INT --> SOLVE["solve_ross_macdonald_multi_site()"]
        RM["ross_macdonald_ode()"] -.-> SOLVE
        SOLVE --> CMP["compute_mechanistic_prediction()"]
    end

    subgraph S2["STAGE 2: GP + Dual Likelihood"]
        CMP -->|"I*(s,t)"| FIT["fit_epiwave_gp()"]
        CASES["Case counts"] --> FIT
        PREV["Prevalence surveys"] --> FIT
        BK["build_gp_kernel()"] -.-> FIT
        FIT --> MCMC["greta::mcmc()"]
    end

    subgraph VAL["VALIDATION"]
        MCMC --> EPS["extract_posterior_summary()"]
        MCMC --> CPM["compute_performance_metrics()"]
    end

    style S1 fill:#EBF5FB,stroke:#2E75B6,stroke-width:2px
    style S2 fill:#E8F8F5,stroke:#27AE60,stroke-width:2px
    style VAL fill:#F5EEF8,stroke:#8E44AD,stroke-width:2px
```

<div class="mt-1 text-xs grid grid-cols-3 gap-2">

<div class="p-1 bg-blue-50 rounded text-center border border-blue-200">
<strong>Stage 1</strong>: 7 functions<br/>
<code>deSolve</code> + <code>approxfun()</code>
</div>

<div class="p-1 bg-green-50 rounded text-center border border-green-200">
<strong>Stage 2</strong>: 4 functions<br/>
<code>greta</code> + <code>greta.gp</code> / TensorFlow HMC
</div>

<div class="p-1 bg-purple-50 rounded text-center border border-purple-200">
<strong>Validation</strong>: 4 functions<br/>
Sim-estimation + metrics
</div>

</div>

<div class="mt-1 text-xs opacity-70 text-center">

`simulate_and_estimate()` orchestrates the full pipeline: data generation, model comparison (GP+offset vs I\*=0 standard geostatistical), and diagnostic plots

</div>

<!--
This is the full architecture of our R implementation — 15 functions total.

Stage 1 has 7 functions. Three parameter generators — get_fixed_m, get_fixed_a, get_fixed_g — each produce time-by-site matrices from Vector Atlas data or temperature-dependent defaults. apply_interventions adjusts m, a, g for ITN/IRS. solve_ross_macdonald_multi_site solves the ODE per site. compute_mechanistic_prediction computes I* = m*a*b*z (infection incidence rate, no population).

Stage 2 has 5 functions. build_gp_kernel constructs the spatial Matern 5/2 kernel. ar1 applies AR(1) temporal correlation (ported from epiwave.mapping). simulate_gp_residuals and simulate_prevalence_surveys generate synthetic data. fit_epiwave_gp builds the GP model with dual Poisson plus Binomial likelihood. Population enters the Poisson likelihood, not I*.

The validation layer has simulate_and_estimate as the orchestrator, plus extract_posterior_summary and compute_performance_metrics.
-->
