---
layout: default
---

# Why This Framework Wins

<div class="grid grid-cols-2 gap-2 mt-2">

<div class="p-2 rounded-lg text-xs bg-blue-50 border border-blue-200">

### Computational Efficiency
- No MCMC over ODE parameters
- One-time forward ODE solve per site
- Only **2 parameters** for HMC
- **17x speedup** vs joint inference

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
- Mechanistic structural prior
- NegBin residual overdispersion
- Data overrides where needed
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
    VA["Vector Atlas<br/>m, a, g"] --> S1["Stage 1<br/>Ross-Macdonald ODE"]
    TEMP["Temperature<br/>Models"] --> S1
    ITN["Interventions<br/>ITN / IRS"] --> S1
    S1 -->|"I*(s,t)"| S2["Stage 2<br/>NegBin Calibration"]
    CASES["DHIS2<br/>Case Data"] --> S2
    S2 -->|"MCMC<br/>(2 params)"| POST["Posterior<br/>Predictions"]
```

<!--
This summarises the four key advantages of the framework. First, computational efficiency — by fixing entomological parameters and solving ODEs once, we reduce the MCMC problem from 20+ parameters to just 2. Second, biological realism — we're using real entomological data and mechanistic transmission dynamics, not just statistical surfaces. Third, statistical flexibility — the Negative Binomial allows the data to override mechanistic predictions where they're wrong. Fourth, operational feasibility — the modular design means you can update vector data independently from case data, and the framework scales to national mapping.

The diagram at the bottom shows the full data flow — Vector Atlas, temperature models, and intervention data feed into Stage 1, which produces the mechanistic prediction I-star. This flows into Stage 2 along with DHIS2 case data, and MCMC inference on just 2 parameters produces posterior predictions.
-->
