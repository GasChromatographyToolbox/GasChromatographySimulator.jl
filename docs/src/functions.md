# Functions 

```@docs
GasChromatographySimulator.Column
```

```@docs
GasChromatographySimulator.Column(L, d, df, sp, gas)
```

```@docs
GasChromatographySimulator.Program
```

```@docs
GasChromatographySimulator.Program(time_steps::Array{<:Real, 1}, temp_steps::Array{<:Real, 1}, pin_steps::Array{<:Real, 1}, pout_steps::Array{<:Real, 1}, L)
```

```@docs
GasChromatographySimulator.Program(time_steps::Array{<:Real, 1}, temp_steps::Array{<:Real, 1}, pin_steps::Array{<:Real, 1}, pout_steps::Array{<:Real, 1}, a_gf::Array{<:Real, 2}, Tcontrol, L)
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

## Helper

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

## Physical Model
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

## Solving
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

## Results
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

```@docs
GasChromatographySimulator.plot_chromatogram
```

```@docs
GasChromatographySimulator.plot_chromatogram!
```

```@docs
GasChromatographySimulator.plot_flow
```

```@docs
GasChromatographySimulator.plot_pressure
```

```@docs
GasChromatographySimulator.plot_temperature
```

```@index
```