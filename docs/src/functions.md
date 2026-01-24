# Functions 

## Structures and Constructors

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
GasChromatographySimulator.solve_separate_multithreads
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

## Plot functions

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

```@docs
GasChromatographySimulator.local_plots
```

```@docs
GasChromatographySimulator.velocity
```

## UI functions for notebooks

```@docs
GasChromatographySimulator.UI_Column
```

```@docs
GasChromatographySimulator.UI_Program
```

```@docs
GasChromatographySimulator.UI_Substance
```

```@docs
GasChromatographySimulator.UI_Options
```

```@docs
GasChromatographySimulator.setting_prog
```

## Helper functions

```@docs
GasChromatographySimulator.temperature_interpolation
```

```@docs
GasChromatographySimulator.steps_interpolation
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

```@docs
GasChromatographySimulator.common
```

```@docs
GasChromatographySimulator.compare_peaklist
```

```@docs
GasChromatographySimulator.compare_measurement_simulation
```

```@docs
GasChromatographySimulator.conventional_program(CP; time_unit="min")
```

```@docs
GasChromatographySimulator.temperature_program(time_steps, value_steps; time_unit="min")
```

```@docs
GasChromatographySimulator.common_time_steps
```

```@docs
GasChromatographySimulator.new_value_steps
```

```@index
```