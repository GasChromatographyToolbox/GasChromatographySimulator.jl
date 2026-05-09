# Changelog

All notable changes to this project will be documented in this file.

## Unreleased

### Added
- Validation helpers for canonical string normalization and numeric sanity checks across constructors (positive/non-negative and finite value guards).
- Finite-value checks in `Substance` construction for core thermodynamic and timing parameters.

### Changed
- `Options` now performs case-insensitive normalization for `Tcontrol`, `vis`, and `control`, validates positive tolerances (`abstol`, `reltol`, `k_th`), and applies a light solver-object check for `alg` (with warning for non-recommended algorithms).
- `Column` constructors now validate `L > 0`, `d > 0`, `df >= 0` and normalize `gas` case-insensitively (`He`, `H2`, `N2`).
- `Program` constructors now validate step-array finiteness/non-negativity, `a_gf` shape, positive `L`, normalized `Tcontrol`, normalized `time_unit` (`min` or `s`), and stricter `pout` handling (`vacuum` / `atmosphere` / non-negative numeric).
- `Parameters` now emits a warning (instead of throwing an error) when duplicate substance names are provided in `sub`, while still rejecting empty `sub`.

### Tests
- Added regression tests for `Options` normalization/validation and warning behavior for non-recommended algorithms.
- Added regression tests for `Column`, `Program`, `Parameters`, and `Substance` constructor validation paths, including duplicate-name warning assertions.

## 0.6.0 - 2026-01-24

### Breaking
- Uncertainty propagation via `Measurements.jl` is currently non-functional (known issue since 0.5.7/0.5.8).

### Added
- Benchmarking utilities for ODE solver comparisons, including `scripts/compare_ode_algorithms.jl`, a new docs page, and benchmark assets.
- Reexported additional OrdinaryDiffEq solvers: `Tsit5`, `Vern9`, `BS5`, `DP5`.
- `BenchmarkTools` dependency for ODE solver benchmarking scripts.
- Optional `label` argument for `plot_chromatogram` and `plot_chromatogram!` to support legends.

### Changed
- ODE solver recommendations and tolerance guidance in `Options` docstrings.
- `compare_measurement_simulation` now accepts measured retention time columns named `:RT` or `:tR`.
- Documentation updates across usage, functions, examples, references, and additional pages, plus new ODE solver guidance.

### Fixed
- `solve_separate` uses the correct length of `Tchar`.
- Reduced redundant residency computations in `odesystem_r!` and `peakode` for better performance.
- `peakode` now accepts an optional precomputed residency value for reuse.
