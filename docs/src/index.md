# GasChromatographySimulator.jl Documentation

```@contents
```

## Introduction
_Aim of this package_

## Installation

## Example

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
