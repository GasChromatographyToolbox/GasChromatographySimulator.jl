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
using DataFrames, CSV # hide
using Plots # hide
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
db_dataframe = DataFrame(CSV.File("../../data/Database_Leppert2020b.csv", header=1, silencewarnings=true))
```
and extract all the names of the substances with:
```@example ex_meas
solutes = GasChromatographySimulator.all_solutes(col.sp, db_dataframe)
```
The data for all solutes is finally loaded with:
```@example ex_meas 
t₀ = zeros(length(solutes))
τ₀ = zeros(length(solutes))
sub = GasChromatographySimulator.load_solute_database(db_dataframe, col.sp, col.gas, solutes, t₀, τ₀)
```
The injection is assumed to be ideal with initial peak widths `τ₀` of 0 seconds and occuring at the beginning of the temperature program (`t₀` of 0 seconds).

The parameters are combined:
```@example ex_meas
par = GasChromatographySimulator.Parameters(col, prog_D, sub, opt)
```

The temperature program and the pressure/flow program can be plotted:
```@example ex_meas
p_flow = GasChromatographySimulator.plot_flow(par)
p_press = GasChromatographySimulator.plot_pressure(par)
p_temp = GasChromatographySimulator.plot_temperature(par)
l = @layout([a{0.65w} [b; c]])
p_TpF = plot(p_temp, p_press, p_flow, layout=l)
```

And the simulation is run:
```@example ex_meas
peaklist, sol = GasChromatographySimulator.simulate(par)
peaklist
```

The file [`Leppert2020b_measured_RT_progD.csv`](https://github.com/JanLeppert/GasChromatographySimulator.jl/blob/main/data/measurements/Leppert2020b_measured_RT_progD.csv) contains the retention times and peak widths (as standard deviations) from the measured chromatogram.
```@example ex_meas
measurement_D = DataFrame(CSV.File("../../data/measurements/Leppert2020b_measured_RT_progD.csv", header=1, silencewarnings=true))
measurement_D[!, 2] = measurement_D[!, 2] .* 60.0 # conversion from min -> s
rename!(measurement_D, [:Name, :tR, :τR])
```

The simulated and measured separations can be compared by comparing the peak lists:
```@example ex_meas
compare = GasChromatographySimulator.compare_peaklist(measurement_D, peaklist)
```
or by comparing the chromatograms:
```@example ex_meas
chrom_D = DataFrame(CSV.File("../../data/measurements/Leppert2020b_measured_Chrom_progD.csv", header=1, silencewarnings=true))
p_chrom, t, chrom = GasChromatographySimulator.plot_chromatogram(peaklist, (0.0, round(chrom_D[end,1];sigdigits=2)); annotation=false, number=true, mirror=true, offset=0.0)
plot!(p_chrom, chrom_D[!,1], chrom_D[!,2].*400.0.+0.1)
ylims!(-1.6,1.6)
xlims!(0.0,round(chrom_D[end,1];sigdigits=2))
p_chrom
```
### Thermal gradient GC

The following example of a thermal gradient GC is the example `medium gradient` from [`[8]`](https://janleppert.github.io/GasChromatographySimulator.jl/dev/references/#References). 

Standard options are used:
```@example ex_meas
opt_tg = GasChromatographySimulator.Options()
```

And the column is defined as:
```@example ex_meas
col_tg = GasChromatographySimulator.Column(2.05, 0.104e-3, 0.104e-6, "FS5ms", "He")
```

The program is taken from the measured temperatures and pressures during the GC run, stored in the file [`Leppert2020b_prog_settings_med_gradient_x90.csv`](https://github.com/JanLeppert/GasChromatographySimulator.jl/blob/main/data/measurements/Leppert2020b_prog_settings_med_gradient_x90.csv):
```@example ex_meas
prog_settings = DataFrame(CSV.File("../../data/measurements/Leppert2020b_prog_settings_med_gradient_x90.csv", header=1, silencewarnings=true))
# use only every 20th measurement
time = cumsum(prog_settings.Deltat)[1:20:end]
time_steps = Array{Float64}(undef, length(time))
for i=2:length(time)
    time_steps[i] = time[i]-time[i-1]
end
time_steps[1] = 0.0
temp_steps = prog_settings.T[1:20:end]
ΔT_steps = prog_settings.DeltaT[1:20:end]
pin_steps = prog_settings.pinj[1:20:end].*1000.0 .+ 101300.0
pout_steps = prog_settings.pdet[1:20:end].*1000.0
α_steps = -3.0.*ones(length(ΔT_steps))
x₀_steps = zeros(length(ΔT_steps))
L₀_steps = col_tg.L.*ones(length(ΔT_steps))
prog_med_grad = GasChromatographySimulator.Program(time_steps, temp_steps, pin_steps, pout_steps, ΔT_steps, x₀_steps, L₀_steps, α_steps, "outlet", col_tg.L)
```

The same solutes are used as in the previous example.

The parameters are combined:
```@example ex_meas
par_tg = GasChromatographySimulator.Parameters(col_tg, prog_med_grad, sub, opt_tg)
```

The temperature program and the pressure/flow program can be plotted:
```@example ex_meas
p_flow_tg = GasChromatographySimulator.plot_flow(par_tg)
p_press_tg = GasChromatographySimulator.plot_pressure(par_tg)
p_temp_tg = GasChromatographySimulator.plot_temperature(par_tg)
l = @layout([a{0.65w} [b; c]])
p_TpF_tg = plot(p_temp_tg, p_press_tg, p_flow_tg, layout=l)
```

And the simulation is run:
```@example ex_meas
peaklist_tg, sol_tg = GasChromatographySimulator.simulate(par_tg)
peaklist_tg
```

The file [`Leppert2020b_measured_RT_med_gradient.csv`](https://github.com/JanLeppert/GasChromatographySimulator.jl/blob/main/data/measurements/Leppert2020b_measured_RT_med_gradient.csv) contains the retention times and peak widths (as standard deviations) from the measured chromatogram.
```@example ex_meas
measurement_tg = DataFrame(CSV.File("../../data/measurements/Leppert2020b_measured_RT_med_gradient.csv", header=1, silencewarnings=true))
measurement_tg[!, 3] = measurement_tg[!, 3] ./ 1000.0 # conversion from ms -> s
rename!(measurement_tg, [:Name, :tR, :τR])
```

The simulated and measured separations can be compared by comparing the peak lists:
```@example ex_meas
compare_tg = GasChromatographySimulator.compare_peaklist(measurement_tg, peaklist_tg)
```
or by comparing the chromatograms:
```@example ex_meas
chrom_tg = DataFrame(CSV.File("../../data/measurements/Leppert2020b_measured_Chrom_med_gradient_x90.csv", header=1, silencewarnings=true))
p_chrom_tg, t_, chrom_ = GasChromatographySimulator.plot_chromatogram(peaklist_tg, (0.0, 55.0); annotation=false, number=true, mirror=true, offset=0.0)
plot!(p_chrom_tg, chrom_tg[!,1], chrom_tg[!,2].*8e-5)
ylims!(-13,13)
xlims!(0.0,55.0)
p_chrom_tg
```