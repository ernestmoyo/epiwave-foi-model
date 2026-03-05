---
layout: default
---

# Code Architecture: `epiwave-foi-model.R`

11 functions in ~580 lines — modular, two-stage pipeline

```mermaid {scale: 0.55}
graph TD
    subgraph S1["STAGE 1: Mechanistic Engine"]
        direction TB
        M["get_fixed_m()<br/><i>mosquito density ratio</i>"] --> INT["apply_interventions()<br/><i>ITN / IRS effects</i>"]
        A["get_fixed_a()<br/><i>biting rate</i>"] --> INT
        G["get_fixed_g()<br/><i>mortality rate</i>"] --> INT
        INT -->|"m, a, g adjusted"| SOLVE["solve_ross_macdonald_multi_site()<br/><i>deSolve ODE per site</i>"]
        RM["ross_macdonald_ode()<br/><i>dx/dt, dz/dt</i>"] -.->|"called internally"| SOLVE
        SOLVE -->|"x(t,s), z(t,s)"| CMP["compute_mechanistic_prediction()<br/><i>I* = m * a * b * z * N</i>"]
    end

    subgraph S2["STAGE 2: Bayesian Calibration"]
        direction TB
        CMP -->|"I*(s,t)"| FIT["fit_epiwave_with_offset()<br/><i>greta NegBin model</i>"]
        CASES["observed_cases<br/><i>DHIS2 / simulated</i>"] --> FIT
        FIT -->|"greta model"| MCMC["greta::mcmc()<br/><i>HMC sampling</i>"]
    end

    subgraph VAL["VALIDATION & METRICS"]
        direction TB
        MCMC --> EPS["extract_posterior_summary()"]
        MCMC --> CPM["compute_performance_metrics()"]
    end

    SIM["simulate_and_estimate()<br/><i>orchestrator: data gen + MCMC + plots</i>"] -.->|"calls all above"| M
    SIM -.-> FIT

    style S1 fill:#EBF5FB,stroke:#2E75B6,stroke-width:2px
    style S2 fill:#E8F8F5,stroke:#27AE60,stroke-width:2px
    style VAL fill:#F5EEF8,stroke:#8E44AD,stroke-width:2px
    style SIM fill:#FEF9E7,stroke:#F39C12,stroke-width:2px
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

<!--
This is the full architecture of our R implementation. The code is 580 lines total — deliberately compact.

Stage 1 has 6 functions. Three parameter generators — get_fixed_m, get_fixed_a, get_fixed_g — each produce time-by-site matrices from either Vector Atlas data or temperature-dependent defaults. These feed into apply_interventions, which adjusts m, a, and g based on ITN and IRS coverage using the Griffin 2010 framework. The adjusted parameters then go to solve_ross_macdonald_multi_site, which calls ross_macdonald_ode internally via deSolve for each site. Finally, compute_mechanistic_prediction combines the ODE output with population data to produce I-star.

Stage 2 has just one core function — fit_epiwave_with_offset — which builds a greta model with the Negative Binomial likelihood and the mechanistic offset. It returns a model object ready for HMC sampling.

The validation layer has simulate_and_estimate as the orchestrator — it generates synthetic data, runs both WITH and WITHOUT models, and produces all 5 diagnostic plots. Plus two helper functions for posterior summaries and performance metrics.

The key design principle is modularity — each function has a single responsibility, and the stages are completely decoupled. You can swap in real Vector Atlas data for Stage 1 without touching Stage 2.
-->
