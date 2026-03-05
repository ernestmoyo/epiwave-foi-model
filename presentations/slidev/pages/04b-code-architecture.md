---
layout: default
---

# Code Architecture: `epiwave-foi-model.R`

11 functions in ~580 lines — modular, two-stage pipeline

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

    subgraph S2["STAGE 2: Bayesian Calibration"]
        CMP -->|"I*(s,t)"| FIT["fit_epiwave_with_offset()"]
        CASES["observed_cases"] --> FIT
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
<strong>Stage 1</strong>: 6 functions<br/>
<code>deSolve</code> + <code>approxfun()</code>
</div>

<div class="p-1 bg-green-50 rounded text-center border border-green-200">
<strong>Stage 2</strong>: 1 function<br/>
<code>greta</code> / TensorFlow HMC
</div>

<div class="p-1 bg-purple-50 rounded text-center border border-purple-200">
<strong>Validation</strong>: 4 functions<br/>
Sim-estimation + metrics
</div>

</div>

<div class="mt-1 text-xs opacity-70 text-center">

`simulate_and_estimate()` orchestrates the full pipeline: data generation, both MCMC models, and all 5 diagnostic plots

</div>

<!--
This is the full architecture of our R implementation. The code is 580 lines total — deliberately compact.

Stage 1 has 6 functions. Three parameter generators — get_fixed_m, get_fixed_a, get_fixed_g — each produce time-by-site matrices from either Vector Atlas data or temperature-dependent defaults. These feed into apply_interventions, which adjusts m, a, and g based on ITN and IRS coverage. The adjusted parameters go to solve_ross_macdonald_multi_site, which calls ross_macdonald_ode internally via deSolve for each site. Finally, compute_mechanistic_prediction combines the ODE output with population data to produce I-star.

Stage 2 has just one core function — fit_epiwave_with_offset — which builds a greta model with the Negative Binomial likelihood and the mechanistic offset.

The validation layer has simulate_and_estimate as the orchestrator, plus two helper functions for posterior summaries and performance metrics.
-->
