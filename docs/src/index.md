# GasChromatographySimulator.jl Documentation



## Introduction

The package GasChromatographySimulator.jl simulates the separation of different substances (solutes) in a [gas chromatographic (GC) system](https://en.wikipedia.org/wiki/Gas_chromatography). The simulation uses ordinary differential equations (ODE) to model the migration ``t(x)`` of a solute through the GC system and the development of the peak variance during this migration ``τ²(x)``.

Beside a temperature program (change of the temperature of the GC-system with time ``T(t)``) and a pressure/flow program (change of inlet ``p_i(t)`` and/or outlet pressure ``p_o(t)`` with time, resp. change of the flow with time ``F(t)``), a thermal gradient (non-uniform change of the temperature along the GC column, ``T(x)``) can be added. Also a non-uniform thickness of the stationary phase ``d_f(x)`` and a non-uniform column diameter ``d(x)`` can be defined. 

The interaction between the substances and the stationary phase of the GC-system is described by a thermodynamic model (K-centric thermodynamic parameters [`[6]`](https://janleppert.github.io/GasChromatographySimulator.jl/dev/references/#References))

![Chromatogram](https://i.ibb.co/HF3gM5r/Chromatogram.png)

For further details see [`[8]`](https://janleppert.github.io/GasChromatographySimulator.jl/dev/references/#References).

The simulation uses the following packages:
- [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl)
- [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl)
- [QuadGK.jl](https://github.com/JuliaMath/QuadGK.jl)
- [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl)
- [Pluto.jl](https://github.com/fonsp/Pluto.jl)

The manual is structured as followed:

```@contents
Pages = [
    "installation.md",
    "usage.md",
    "functions.md",
    "references.md"
    ]
Depth = 2
```


