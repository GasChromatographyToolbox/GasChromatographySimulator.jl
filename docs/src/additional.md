# Additional features

## Differentiability

The solutions of the ODEs, the retention times and the peak width, can be differentiated for different parameters of the system, e.g., the length ``L`` or retention parameters ``T_{\text{char}}``, using automatic differentiation (`ForwardDiff.jl`).

!!! warning "Attention!"

    The functionality of differentiability at the moment is intended only for non-gradient cases (no spatial thermal gradient, no variation of column diameter or film thickness along the column). Therefore, the option `ng` should be set to `ng = true`, see  [`GasChromatographySimulator.Options`](@ref).     

Example code:

Using the solution of the  ODE-system, see  [`GasChromatographySimulator.solving_odesystem_r`](@ref) for a substance, the resulting retention time `tR` and peak variance `τR²` can be differentiated for selected parameters, e.g. the three retention parameters `Tchar`, `θchar` and `ΔCp`. A new function depending on these parameter (`x = [Tchar, θchar, ΔCp]`) is created. 
```@example diff
using GasChromatographySimulator
using ForwardDiff
L = 30.0; d = 0.25e-3; df = 0.25e-6; gas = "He";
prog = GasChromatographySimulator.Program([40.0, 3.0, 10.0, 300.0, 5.0], [300000.0, 3.0, (450000.0-300000.0)/(300.0-40.0)*10.0, 450000.0, 5.0], L);
φ₀ = 1e-3; Cag = 1e-6; t₀ = 0.0; τ₀ = 0.0;
opt = GasChromatographySimulator.Options(ng=true);  
tR_τR2_RP(x) = GasChromatographySimulator.solving_odesystem_r(L, d, df, gas, prog.T_itp, prog.Fpin_itp, prog.pout_itp, x[1], x[2], x[3], φ₀, Cag, t₀, τ₀, opt).u[end];
tR_τR2_RP([400.0, 30.0, 100.0])
```
This function can be differentiated for these parameters using an automatic differentiation, like [`ForwardDiff.jl`](https://github.com/JuliaDiff/ForwardDiff.jl). 
```@example diff
∂_tR_τR2_RP(x) = ForwardDiff.jacobian(tR_τR2_RP, x)
nothing # hide
```
This function calculated the derivative for the selected value:
```@example diff
∂_tR_τR2_RP([400.0, 30.0, 100.0])
```

This functionality can be used for sensitivity analysis and can be helpful for the use of optimization methods.
## Uncertainty

A uncertainty for the solutions of the ODEs can be calculated based on uncertainties of parameters of the system, e.g., column diameter ``d`` or retention parameters ``\Delta C_\text{p}``, using the Julia package `Measurements.jl`.

!!! warning "Attention!"

    The functionality of uncertainties at the moment is intended only for non-gradient cases (no spatial thermal gradient, no variation of column diameter or film thickness along the column). Therefore, the option `ng` should be set to `ng = true`, see  [`GasChromatographySimulator.Options`](@ref).     

Example code

For the example the values of the column dimensions, length `L`, diameter `d` and film thickness `df` are subjected to uncertainties.
```@example uncertainty
using GasChromatographySimulator
using Measurements
col = GasChromatographySimulator.Column(30.0±1.0, (0.25±0.01)*1e-3, (0.25±0.01)*1e-6, "Rxi5SilMS", "He")
prog = GasChromatographySimulator.Program([40.0, 3.0, 10.0, 300.0, 5.0], [300000.0, 3.0, (450000.0-300000.0)/(300.0-40.0)*10.0, 450000.0, 5.0], col.L)
solutes = ["C10", "C11", "C12", "2-Octanol", "2-Octanone"]
t₀ = fill(0.0±0.0, length(solutes))
τ₀ = fill(0.0±0.0, length(solutes))
sub = GasChromatographySimulator.load_solute_database("../../data", "Database_test.csv", 
                                                        col.sp,
                                                        col.gas,
                                                        solutes,
                                                        t₀,
                                                        τ₀)
opt = GasChromatographySimulator.Options(ng=true) 
par = GasChromatographySimulator.Parameters(col, prog, sub, opt)                                                       
nothing # hide
```
The initial time `t₀` and peak width `τ₀` must also be defined with an uncertainty. 
```@example uncertainty
pl, sol = GasChromatographySimulator.simulate(par)
```
The result is expressed also with uncertainties, e.g. the peaklist:
```@example uncertainty
pl # hide
```
The chromatogram can be plotted with markers for the uncertainties:
PLOT_EXAMPLE