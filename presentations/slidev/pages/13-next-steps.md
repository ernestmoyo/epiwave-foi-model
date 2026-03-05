---
layout: default
---

# Next Steps

From simulation validation to production-scale Angola mapping

<div class="grid grid-cols-2 gap-2 mt-2">

<div class="p-2 rounded-lg bg-blue-50 border-l-4 border-blue-500 text-xs">

### 1. Angola DHIS2 Integration
18 provinces × 60+ months of routine case data from NMCP via DHIS2
- Quality checks for reporting completeness
- Outlier detection and data cleaning
- Monthly aggregated case counts

</div>

<div class="p-2 rounded-lg bg-green-50 border-l-4 border-green-500 text-xs">

### 2. Vector Atlas Parameters
Empirical estimates for $m$, $a$, $g$ from entomological surveys
- Map survey-level vector parameters to province-level inputs
- Temperature-dependent parameter models as fallback
- Incorporate insecticide resistance data

</div>

<div class="p-2 rounded-lg bg-purple-50 border-l-4 border-purple-500 text-xs">

### 3. INLA/SPDE Spatial Extension
Re-introduce spatial correlation via Matern SPDE on mesh (R-INLA)
- Scales to **1000+ spatial units** without Cholesky issues
- Replaces failed GP approach with sparse precision matrices
- Proven technology in malaria mapping (MAP, Bhatt et al.)

</div>

<div class="p-2 rounded-lg bg-amber-50 border-l-4 border-amber-500 text-xs">

### 4. INFO Document & Review
Complete framework validation write-up for supervisors
- Simulation results and convergence diagnostics
- Comparison with existing geostatistical approaches
- Manuscript preparation for Objective 2 paper

</div>

</div>

<!--
Looking ahead, there are four key next steps. First, integrating real data from Angola's DHIS2 system — 18 provinces with over 60 months of routine case data. This will be the first real-world test of the framework.

Second, obtaining empirical Vector Atlas parameters for Angola — mapping entomological survey data to province-level inputs for m, a, and g. Where survey data aren't available, we'll use temperature-dependent models as fallback.

Third, the INLA/SPDE extension. Now that we've validated the offset approach works with the Negative Binomial, we can re-introduce spatial correlation using the Matern SPDE on a mesh via R-INLA. This avoids the Cholesky issues that killed the GP approach because SPDE uses sparse precision matrices instead of dense covariance matrices. This is proven technology — it's what MAP uses.

Fourth, completing the framework documentation and beginning manuscript preparation for the Objective 2 paper.
-->
