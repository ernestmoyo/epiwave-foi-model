# Correspondence & Issues Log

This folder tracks issues, decisions, and communications related to the EpiWave FOI Model development.

## Structure

- `issues/` — Technical issues raised, investigated, and resolved
- `decisions/` — Key architectural decisions shared with supervisors
- `updates/` — Code update summaries sent to collaborators

## Issues Index

### Architecture & Design
| File | Date | Topic | Status |
|------|------|-------|--------|
| `2026-01-24_two_stage_architecture.md` | 2026-01-24 | Two-stage architecture (the brainwave) | Resolved |
| `2026-03-05_stage2_gp_to_negbin.md` | 2026-03-05 | GP replaced with NegBin in Stage 2 | Resolved |

### Computational Performance
| File | Date | Topic | Status |
|------|------|-------|--------|
| `2026-01-24_ode_speed_test_r_vs_tf.md` | 2026-01-24 | R/deSolve vs TF DormandPrince speed | Resolved |
| `GH-004_ode_solver_speed_benchmarks.md` | 2026-01-24 | GitHub Issue #4 — ODE benchmarks | Resolved |
| `GH-005_tf_ode_solver_integration.md` | 2026-01-24 | GitHub Issue #5 — TF ODE integration | Resolved |

### Statistical / Convergence
| File | Date | Topic | Status |
|------|------|-------|--------|
| `2026-03-02_alpha_gamma_nonidentifiability.md` | 2026-03-02 | alpha/gamma non-identifiability | Resolved |
| `2026-03-05_ess_convergence_phi_reparameterisation.md` | 2026-03-05 | ESS=82 with Beta(phi), fixed with log_size | Resolved |

### Environment / Setup
| File | Date | Topic | Status |
|------|------|-------|--------|
| `2025-09-18_greta_offline_bug.md` | 2025-09-18 | Greta fails when offline (UV_OFFLINE fix) | Resolved |
| `2026-02-25_pkg_resources_python_environment.md` | 2026-02-25 | pkg_resources missing, self-healing loader | Resolved |

## Reading Order (for understanding the journey)

1. **ODE speed test** — discovered TF is too slow for MCMC
2. **Two-stage architecture** — Nick's brainwave to fix the speed problem
3. **GP to NegBin** — why Nick's original GP Stage 2 didn't work computationally
4. **alpha/gamma non-identifiability** — why we merged to log_rate
5. **ESS convergence** — why we reparameterised from phi to log_size
6. Environment issues (offline bug, pkg_resources) — practical blockers resolved

## How to use

When sharing updated code or raising a point with supervisors, create a dated markdown file documenting:
1. What changed and why
2. Evidence (experiment results, diagnostics)
3. References to code/experiments
