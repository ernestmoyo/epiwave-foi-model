---
layout: default
---

# Next Steps

From simulation validation to production-scale Angola mapping

<div class="grid grid-cols-2 gap-4 mt-4">

<div class="p-3 rounded-lg bg-blue-50 border-l-4 border-blue-500 text-sm">

### 1. Angola DHIS2 Integration
18 provinces × 60+ months of routine case data from NMCP via DHIS2
- Quality checks for reporting completeness
- Outlier detection and data cleaning
- Monthly aggregated case counts

</div>

<div class="p-3 rounded-lg bg-green-50 border-l-4 border-green-500 text-sm">

### 2. Vector Atlas Parameters
Empirical estimates for $m$, $a$, $g$ from entomological surveys
- Map survey-level vector parameters to province-level inputs
- Temperature-dependent parameter models as fallback
- Incorporate insecticide resistance data

</div>

</div>

<!--
Looking ahead, there are four key next steps. First, integrating real data from Angola's DHIS2 system — 18 provinces with over 60 months of routine case data. This will be the first real-world test of the framework.

Second, obtaining empirical Vector Atlas parameters for Angola — mapping entomological survey data to province-level inputs for m, a, and g. Where survey data aren't available, we'll use temperature-dependent models as fallback.

Third, the INLA/SPDE extension. Now that we've validated the offset approach works with the Negative Binomial, we can re-introduce spatial correlation using the Matern SPDE on a mesh via R-INLA. This avoids the Cholesky issues that killed the GP approach because SPDE uses sparse precision matrices instead of dense covariance matrices. This is proven technology — it's what MAP uses.

Fourth, completing the framework documentation and beginning manuscript preparation for the Objective 2 paper.
-->
