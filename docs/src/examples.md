# Examples

## Notebooks

In the folder [notebooks](https://github.com/JanLeppert/GasChromatographySimulator.jl/tree/main/notebooks) several notebooks, using [Pluto.jl](https://github.com/fonsp/Pluto.jl), for the simulation of GC-systems are available. 

To use these notebooks [Julia, v1.6 or above,](https://julialang.org/downloads/#current_stable_release) must be installed and **Pluto** must be added:

```julia
julia> ]
(v1.7) pkg> add Pluto
```

To run Pluto, use the following commands:

```julia
julia> using Pluto
julia> Pluto.run()
```

Pluto will open your browser. In the field `Open from file` the URL of a notebook or the path to a locally downloaded notebook can be insert and the notebook will open and load the necessary packages. 

### Overview of notebooks

- `simulation_conventional_GC.jl` - simulation of a conventional GC-system (constant temperature, constant diameter and constant film thickness along the column) and outlet pressure as "vacuum" or "atmospheric", now with option of flow or pressure control and temperature program notation in the typical form used in commercial GC software (temperature levels, holding times and heating ramps) 
- `simulation_conventional_GC_load_2dbs.jl` - simulation of a conventional GC-system (constant temperature, constant diameter and constant film thickness along the column) and loading of up to two different substance databases and simulation of the common substances with the same GC-system and comparing the result. Also, an option is given, to load measured retention times and compare these to the simulations. Same setting for programs as in `simulation_conventional_GC.jl`
- `simulation_example.jl` - general example of simulation of a GC-system with optional thermal gradient (exponential/linear model of temperature change along the column) and constant diameter and constant film thickness along the column. 
- `simulation_example_input_gradient_function.jl` - simulation of a GC-system with optional thermal gradient where the temperature change along the column is defined by a user-defined equation (cosine-function as example)

## Simulation of measurements

Two measurements from [`[8]`](https://janleppert.github.io/GasChromatographySimulator.jl/dev/references/#References) will be simulated and compared to the measured chromatograms. The n-alkanes from n-nonane (C9) to n-triacontane (C30) are separated in a conventional GC and a thermal gradient GC . The database with the thermodynamic parameters is [`Database_Leppert2020b.csv`](https://github.com/JanLeppert/GasChromatographySimulator.jl/blob/main/data/Database_Leppert2020b.csv).

### Conventional GC

The conventional GC program simulated here is `Prog. D` from [`[8]`](https://janleppert.github.io/GasChromatographySimulator.jl/dev/references/#References), a temperature program with two heating ramps, constant inlet pressure and a flame ionization detector (FID, atmospheric outlet pressure). 

The standard options are used. Only the option `ng` (non-gradient) is changed to `true`. Because the conventional GC does not use non-uniform temperature, diameter or film thickness, the model can be simplified and the calculation of the separation is faster.
```@example ex_meas
using GasChromatographySimulator # hide
opt = GasChromatographySimulator.Options(ng=true)
```

The column is defined as:
```@example ex_meas
col = GasChromatographySimulator.Column(11.18, 0.104e-3, 0.104e-6, "FS5ms", "H2")
```

The program is defined as:
```@example ex_meas
prog_D = GasChromatographySimulator.Program([0.0, 60.0, 1680.0, 60.0, 360.0, 60.0], [40.0, 40.0, 180.0, 180.0, 300.0, 300.0], 411564.0*ones(6), 101300.0.*ones(6), col.L)
```

We want to use all solutes for the stationary phase FS5ms, which are in the database. We load the database into a dataframe:
```@example ex_meas
db_dataframe = DataFrame(CSV.File("/../../data/Database_Leppert2020b.csv", header=1, silencewarnings=true))
```
and extract all the names of the substances with:
```@example ex_meas
solutes = GasChromatographySimulator.all_solutes(col.sp, db_dataframe)
```
The data for all solutes is finally loaded with:
```@example ex_meas 
t₀ = 0.8.*ones(length(solutes))
τ₀ = zeros(length(solutes))
sub = GasChromatographySimulator.load_solute_database(db_dataframe, col.sp, col.gas, solutes, t₀, τ₀)
```
Here the injection time `t₀` of all substances is set to 0.8s, because the autosampler used with the GC starts the temperature program shortly before injection. Otherwise the injection is assumed to be ideal with initial peak widths `τ₀` of 0 seconds.

The parameters are combined:
```@example ex_meas
par = GasChromatographySimulator.Parameters(col, prog_D, sub, opt)
```
And the simulation is run:
```@example ex_meas
peaklist, sol = GasChromatographySimulator.simulate(par)
```

THe file [`Leppert2020b_measured_RT_progD.csv`](https://github.com/JanLeppert/GasChromatographySimulator.jl/blob/main/data/Leppert2020b_measured_RT_progD.csv) contains the retention times and peak widths (as standard deviations) from the measured chromatogram.
```@example ex_meas
measurement_D = DataFrame(CSV.File("/../../data/Leppert2020b_measured_RT_progD.csv", header=1, silencewarnings=true))
measurement_D = measurement_D[!, 2] .* 60.0 # conversion from min -> s
rename!(measurement_D, [:Name, :tR, :τR])
```

The simulated and measured separations can be compared by comparing the peak lists:
```@example ex_meas
compare = GasChromatographySimulator.compare_peaklist(measurement_D, peaklist)
```
or by comparing the chromatograms:
**PlotChromatograms**

### Thermal gradient GC

**here an example of a simulation of one of the thermal gradient measurements in Leppert2020b** 