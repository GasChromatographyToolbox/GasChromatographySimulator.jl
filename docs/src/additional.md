# Additional features

## Differentiability

The solutions of the ODEs, the retention times and the peak width, can be differentiated for different parameters of the system, e.g., the length ``L`` or retention parameters ``T_{\text{char}}``, using automatic differentiation (`ForwardDiff.jl`).

...

## Uncertainty

A uncertainty for the solutions of the ODEs can be calculated based on uncertainties of parameters of the system, e.g., column diameter ``d`` or retention parameters ``\Delta C_\text{p}``, using the Julia package `Measurements.jl`.

...