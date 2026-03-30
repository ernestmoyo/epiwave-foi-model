---
layout: default
---

# Why This Framework Wins

<div class="grid grid-cols-2 gap-2 mt-2">

<div class="p-2 rounded-lg text-xs bg-blue-50 border border-blue-200">

### Computational Efficiency
- No MCMC over ODE parameters
- One-time forward ODE solve per site
- Only **5 parameters** for HMC
- Sparse GP with inducing points

</div>

<div class="p-2 rounded-lg text-xs bg-green-50 border border-green-200">

### Biological Realism
- Fixed Vector Atlas parameters
- Ross-Macdonald dynamics
- ITN/IRS intervention effects
- Seasonal vector abundance

</div>

<div class="p-2 rounded-lg text-xs bg-purple-50 border border-purple-200">

### Statistical Flexibility
- GP captures spatial residuals
- Dual likelihood (cases + prevalence)
- Data overrides mechanistic model where needed
- Bayesian uncertainty quantification

</div>

<div class="p-2 rounded-lg text-xs bg-amber-50 border border-amber-200">

### Operational Feasibility
- Scales to national mapping
- Modular independent updates
- Interpretable bio parameters
- Counterfactual scenarios

</div>

</div>

```mermaid {scale: 0.4}
graph LR
    VA["Vector Atlas<br/>get_fixed_m/a/g()"] --> INT["apply_interventions()"]
    INT --> S1["solve_ross_macdonald<br/>_multi_site()"]
    S1 -->|"compute_mechanistic<br/>_prediction()"| ISTAR["I*(s,t)"]
    ISTAR --> S2["fit_epiwave_gp()"]
    CASES["DHIS2<br/>Case Data"] --> S2
    PREV["Prevalence<br/>Surveys"] --> S2
    S2 -->|"greta::mcmc()<br/>(5 params)"| POST["Posterior<br/>Predictions"]
```

<!--
This summarises the four key advantages of the framework. First, computational efficiency — by fixing entomological parameters and solving ODEs once, we avoid MCMC over ODE parameters entirely. The sparse GP with inducing points keeps Stage 2 tractable. Second, biological realism — we're using real entomological data and mechanistic transmission dynamics, not just statistical surfaces. Third, statistical flexibility — the GP captures spatially-correlated departures from the mechanistic prediction, and the dual likelihood makes all parameters identifiable. Fourth, operational feasibility — the modular design means you can update vector data independently from case data, and the framework scales to national mapping.

The diagram shows the full data flow — Vector Atlas and intervention data feed into Stage 1, producing I-star. This flows into Stage 2 along with DHIS2 case data and prevalence surveys, and MCMC inference on 5 parameters produces posterior predictions with spatial uncertainty.
-->
