# Issue: Greta Session Fails When Offline

**Date:** 2025-09-18
**Author:** Ernest Moyo (fix by Nick Golding)
**Status:** Resolved — workaround in `R/experiments/tf_ode_test.R`
**Related commit:** `e390b96`

---

## Summary

When working without internet access (common in Dar es Salaam), greta's Python/TensorFlow initialization would fail because the reticulate UV package manager attempted network requests during startup.

## Symptoms

- `library(greta)` hangs or errors when offline
- Error messages related to UV package resolution or network timeouts
- Blocked all development work when internet was unavailable

## Workaround

Nick identified the fix — set the UV_OFFLINE environment variable before loading TensorFlow:

```r
Sys.setenv(UV_OFFLINE = 1)
```

This tells the UV package manager to use only cached packages and skip network requests.

## Where Applied

- `R/experiments/tf_ode_test.R` line 32: `Sys.setenv(UV_OFFLINE = 1)`
- Should be set in `.Rprofile` or before any `library(greta)` call when working offline

## Context

This was discovered during the September 2025 Dar es Salaam intensive (Nick, Gerry, Dave visiting). Working at NM-AIST where internet connectivity is intermittent, this was a practical blocker that needed an immediate fix.

## References

- Commit: `e390b96` — "work around offline greta bug"
- Applied in: `R/experiments/tf_ode_test.R`
