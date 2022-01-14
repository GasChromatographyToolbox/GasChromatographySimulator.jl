# GasChromatographySimulator.jl Documentation

```@contents
```

## Introduction
_Aim of this package_

The package GasChromatographySimulator.jl simulates the separation of different substances (solutes) in a gas chromatographic (GC) system. THe simulation uses ordinary differential equations (ODE) to model the migration (t(z)) of a solute through the GC system and the development of the peak variance during this migration (τ²(z)).

...

For further details see [Leppert.2020b].

This package uses (and reexports) DifferentialEquations.jl [DifferentialEquations.jl] and Interpolations.jl [Interpolations.jl]. Also QuadGK.jl [QuadGK.jl] and ForwardDiff.jl [ForwardDiff.jl] are used.


## Installation

To use GasChromatographicSimulator, you need to install Julia 1.6 or greater first (official Julia website) and than add the package:

```julia
julia> ] add https://github.com/JanLeppert/GasChromatographySimulator.jl
```

## Example

First a simple example.

Advanced example with a database file.

Advanced example with a thermal gradient.

(Advanced example with a user defined gradient function (or diameter function, or film thickness function))

## Pluto Notebook

A Pluto notebook with a simple example is available. [Pluto.jl]

## Functions 

### Structures

```@docs
GasChromatographySimulator.System
```

```@docs
GasChromatographySimulator.Program
```

```@docs
GasChromatographySimulator.Substance
```

```@docs
GasChromatographySimulator.Options
```

```@docs
GasChromatographySimulator.Parameters
```

### Helper

```@docs
GasChromatographySimulator.temperature_interpolation
```

```@docs
GasChromatographySimulator.pressure_interpolation
```

```@docs
GasChromatographySimulator.load_solute_database
```

```@docs
GasChromatographySimulator.all_solutes
```

```@docs
GasChromatographySimulator.diffusivity
```

### Physical Model
```@docs
GasChromatographySimulator.pressure
```

```@docs
GasChromatographySimulator.flow_restriction
```

```@docs
GasChromatographySimulator.viscosity
```

```@docs
GasChromatographySimulator.holdup_time
```

```@docs
GasChromatographySimulator.flow
```

```@docs
GasChromatographySimulator.mobile_phase_residency
```

```@docs
GasChromatographySimulator.residency
```

```@docs
GasChromatographySimulator.retention_factor
```

```@docs
GasChromatographySimulator.plate_height
```

```@docs
GasChromatographySimulator.diffusion_mobile
```

### Solving
```@docs
GasChromatographySimulator.simulate
```

```@docs
GasChromatographySimulator.solve_system_multithreads
```

```@docs
GasChromatographySimulator.solve_multithreads
```

```@docs
GasChromatographySimulator.solving_migration
```

```@docs
GasChromatographySimulator.solving_peakvariance
```

```@docs
GasChromatographySimulator.solving_odesystem_r
```

```@docs
GasChromatographySimulator.odesystem_r!
```

```@docs
GasChromatographySimulator.peakode
```

### Results
```@docs
GasChromatographySimulator.peaklist(sol, par)

GasChromatographySimulator.peaklist(sol, peak, par)
```

```@docs
GasChromatographySimulator.sol_extraction(sol, par)
```

```@docs
GasChromatographySimulator.sol_extraction(sol, peak, par)
```

## Index

```@index
```

## References

[Leppert.2020a]
[DifferentialEquations.jl]
[Interpolations.jl]
[QuadGK.jl]
[ForwardDiff.jl]
[Pluto.jl]