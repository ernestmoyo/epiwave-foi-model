# Issue: pkg_resources Missing in Reticulate Python Environment

**Date:** 2026-02-25
**Author:** Ernest Moyo
**Status:** Resolved — self-healing loader created
**Related code:** `R/greta_setup.R`, `R/pkg_resources_shim.py`

---

## Summary

Greta depends on TensorFlow Probability, which imports `pkg_resources` (part of setuptools). The reticulate-managed Python virtual environment sometimes lacks setuptools, causing greta to fail at import time with:

```
ModuleNotFoundError: No module named 'pkg_resources'
```

This was particularly frustrating because it would appear intermittently — after R restarts, package updates, or on different machines.

## Root Cause

Reticulate's UV-based virtual environment manager creates minimal venvs that may not include setuptools. TensorFlow Probability's import chain requires `pkg_resources.require()` for version checking. Without setuptools installed, the entire greta stack fails.

## Failed Approaches

1. **Manual pip install:** Running `pip install setuptools` in the venv works temporarily but is fragile — any reticulate update or venv recreation breaks it again.

2. **Adding to .Rprofile:** Setting up the environment at startup helps but doesn't work for background jobs or sourced scripts.

## Solution: Self-Healing Loader

Created `R/greta_setup.R` — a robust loader that:

1. Checks if `pkg_resources` is already importable
2. If not, locates `R/pkg_resources_shim.py` (a minimal shim that provides the `require()` function TFP needs)
3. Deploys the shim into the venv's site-packages directory
4. Verifies the fix worked
5. Then loads greta and greta.gp

The shim (`R/pkg_resources_shim.py`) is a minimal Python module that provides just enough of the `pkg_resources` API to satisfy TensorFlow Probability's import requirements.

### Key Design Decisions

- **Self-healing:** Repairs itself automatically on every load, no manual intervention needed
- **Works everywhere:** Interactive sessions, sourced scripts, background jobs, .Rprofile
- **Robust path resolution:** Searches multiple locations for the shim file
- **Non-destructive:** Only deploys if `pkg_resources` is actually missing

## Usage

Add to the top of any script that uses greta:
```r
if (file.exists('R/greta_setup.R')) source('R/greta_setup.R')
```

Or source it from `.Rprofile` for automatic loading.

## References

- Self-healing loader: `R/greta_setup.R`
- Python shim: `R/pkg_resources_shim.py`
- Test framework uses it: `test_framework.R` line 22
