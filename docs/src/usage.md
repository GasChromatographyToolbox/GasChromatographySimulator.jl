# Usage

## Parameters defining a simulation

A GC-system for the simulation is defined by four sets of parameters:

### GC-system

The GC-system parameters `GasChromatographySimulator.System` defines the dimensions of the GC column, length ``L``, diameter ``d`` and film thickness of the stationary phase ``d_f``, all measured in meters, the name of the stationary phase and the name of the mobile phase (with the allowed values "He", "H2" and "N2").

![GC-column](https://i.ibb.co/0y2zqdG/Column.png)

```julia
sys = GasChromatographySimulator.System(4.0, 0.1e-3, 0.1e-6, "Rxi17SilMS", "He")
```

### Program

The program parameters `GasChromatographySimulator.Program` defines the temperature and pressure program for the GC separation.

The definition of the program parameters will be explained by two examples.

#### Without thermal gradient

Without a thermal gradient the temperature is the same at every column position at the same time. This is the normal case for conventional GC. One example of such a program can be achieved by the following method, which constructs the Program-structure:

```julia
prog = GasChromatographySimulator.Program(  [0.0, 60.0, 600.0, 120.0],
                                            [40.0, 40.0, 300.0, 300.0],
                                            [18.0, 18.0, 98.0, 98.0].*1000.0 .+ 101300.0,
                                            [0.0, 0.0, 0.0, 0.0],
                                            sys.L)
```

The first array defines the time steps (in s), the second array defines the temperatures (in °C) at these time steps, the third and fourth array define the inlet and outlet pressures (both in Pa(absolute)) at the time steps. The values of temperature and pressures change linearly between the values defined at the time steps. The following picture shows the resulting temperature and pressure program.

![Program without thermal gradient](https://i.ibb.co/wLdNzzm/T-of-t-and-p-of-t-ng.png)

The first time step is always zero (t₁ = 0.0 s). The following time steps define the time that passes until the next step. In the example the second time step is t₂ = 60 seconds long and in this time the temperature stays constant at 40°C (it changes linearly from T₁ = 40°C to T₂ = 40°C). With the next time step (t₃ = 600 s) the temperature changes from T₂ = 40°C linearly to T₃ = 300°C. In the last time step (t₄ = 120 s) the temperature is again kept constant at 300°C. The pressure program is defined in the same way.

The four arrays for time steps, temperatures and the two pressures must have the same number of elements, otherwise the construction of the Program-structure gives an error message.

#### Thermal gradient
...

...pic_thermal_gradient

### Substances

#### Database

### Additional Options

## Run the simulation

## Evaluate the results
