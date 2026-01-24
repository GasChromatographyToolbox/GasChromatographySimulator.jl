# ODE Solvers

This package solves the GC migration and peak-width ODEs through
`OrdinaryDiffEq.jl`. The solver is selected via `Options(alg=...)`. Common
solvers are re-exported from `GasChromatographySimulator`, so they can be used
without importing `OrdinaryDiffEq` directly.

## Recommended choices

- **`OwrenZen5()` (default)**: Good accuracy and stability for GC simulations.
- **`Tsit5()`**: Modern, efficient 5th-order method. Often faster than
  `OwrenZen5()` while keeping similar accuracy, but validate for your programs.
- **`Vern9()`**: High-accuracy solver. Use for precision-critical studies and
  tighten tolerances.
- **`BS5()`**: Alternative 5th-order method with robust error control.
- **`DP5()`**: Classic 5th-order method. Can be fast, but validate carefully.
- **`OwrenZen3()` / `OwrenZen4()`**: Legacy or speed-focused options.

## How to select a solver

Use `Options(alg=...)` when you create parameters:

```@example ex
using GasChromatographySimulator

opt = GasChromatographySimulator.Options(alg=Tsit5(), abstol=1e-6, reltol=1e-3)
```

For high-accuracy solvers (e.g. `Vern9()`), tighten tolerances:

```@example ex
opt_hi = GasChromatographySimulator.Options(alg=Vern9(), abstol=1e-8, reltol=1e-5)
```

## Validation and benchmarking

The script `scripts/compare_ode_algorithms.jl` compares multiple solvers on both a
long-column case and a short-column gradient case. Use it whenever you change
solver settings or want to confirm numerical equivalence across solvers.

Key guidance:

- Validate solver changes with at least one non-gradient and one gradient case.
- Compare both retention times (`tR`) and peak widths (`τR`).
- If differences exceed your acceptable threshold, keep `OwrenZen5()` or tighten
  tolerances.

### Example benchmark results (local tests)

Results from a representative comparison (17 solutes, "Wax" stationary phase):

- **No gradient (30 m column):**
  - All solvers succeeded.
  - Max differences vs `OwrenZen5()` were small for `Tsit5()`, `Vern9()`,
    `BS5()`, and `OwrenZen4()`.
  - `DP5()` showed larger deviations for some solutes and should be validated
    carefully before use.

- **With thermal gradient (2 m column, stronger ΔT):**
  - All solvers succeeded.
  - Differences vs `OwrenZen5()` were generally smaller than in the no-gradient
    case for most solvers.
  - `Tsit5()` and `Vern9()` showed good agreement with the baseline.

These results are sensitive to the temperature program, column dimensions, and
solute set. Always confirm solver agreement for your specific configuration.

### Benchmark plots

#### Runtime comparison (no gradient)

![Runtime (no gradient)](assets/ode_solver_runtime_no_gradient.png)

#### Max retention time differences (gradient)

![Max |ΔtR| vs baseline (gradient)](assets/ode_solver_maxdiff_gradient.png)

#### Chromatogram overlays (no gradient)

![Chromatogram comparison (no gradient)](assets/ode_solver_chrom_no_gradient.png)

#### Chromatogram overlays (gradient)

![Chromatogram comparison (gradient)](assets/ode_solver_chrom_gradient.png)

## Notes on gradients

Thermal gradients can amplify numerical differences. If you use gradients or
short columns, re-check solver agreement with the gradient test case in
`scripts/compare_ode_algorithms.jl` and tighten tolerances if needed.
