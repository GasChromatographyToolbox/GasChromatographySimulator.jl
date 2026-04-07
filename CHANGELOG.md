# Changelog

All notable changes to this project will be documented in this file.

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
