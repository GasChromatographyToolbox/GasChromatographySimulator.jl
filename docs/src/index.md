# GasChromatographySimulator.jl Documentation



## Introduction

The package GasChromatographySimulator.jl simulates the separation of different substances (solutes) in a [gas chromatographic (GC) system](https://en.wikipedia.org/wiki/Gas_chromatography). The simulation uses ordinary differential equations (ODE) to model the migration ``t(z)`` of a solute through the GC system and the development of the peak variance during this migration ``τ²(z)``.

Beside a temperature program (change of the temperature of the GC-system with time) and a pressure program (change of inlet and/or outlet pressure with time), a thermal gradient (change of the temperature along the GC column) can be added. 

The interaction between the substances and the stationary phase of the GC-system is described by a thermodynamic model (K-centric thermodynamic parameters [Blumberg.2017])

![Chromatogram](https://i.ibb.co/HF3gM5r/Chromatogram.png)

For further details see [Leppert.2020a].

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


