# GasChromatographySimulator.jl Documentation



## Introduction

The package GasChromatographySimulator.jl simulates the separation of different substances (solutes) in a [gas chromatographic (GC) system](https://en.wikipedia.org/wiki/Gas_chromatography). The simulation uses ordinary differential equations (ODE) to model the migration ``t(x)`` of a solute through the GC system and the development of the peak variance during this migration ``τ²(x)``.

Beside a temperature program (change of the temperature of the GC-system with time ``T(t)``) and a pressure/flow program (change of inlet ``p_i(t)`` and/or outlet pressure ``p_o(t)`` with time, resp. change of the flow with time ``F(t)``), a thermal gradient (non-uniform change of the temperature along the GC column, ``T(x)``) can be added. Also a non-uniform thickness of the stationary phase ``d_f(x)`` and a non-uniform column diameter ``d(x)`` can be defined. 

The interaction between the substances and the stationary phase of the GC-system is described by a thermodynamic model (K-centric thermodynamic parameters [`[6]`](https://janleppert.github.io/GasChromatographySimulator.jl/dev/references/#References))

![Chromatogram](https://i.ibb.co/HF3gM5r/Chromatogram.png)

For further details see [`[8]`](https://janleppert.github.io/GasChromatographySimulator.jl/dev/references/#References).

The simulation uses the following packages:
- [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) part of [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl)
- [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl)
- [Integrals.jl](https://github.com/SciML/Integrals.jl)
- [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl)
- [ChemicalIdentifiers.jl](https://github.com/longemen3000/ChemicalIdentifiers.jl)
- [Pluto.jl](https://github.com/fonsp/Pluto.jl)

## Contribution

Please open an issue if you:
- want to report a bug 
- have problems using the package (please first look at the documentation)
- have ideas for new features or ways to improve the usage of this package 

You can contribute (e.g. fix bugs, add new features, add to the documentation) to this package by Pull Request: 
- first discuss your contributions in a new issue
- ensure that all tests pass locally before starting the pull request
- new features should be included in `runtests.jl`
- add description to the pull request, link to corresponding issues by `#` and issue number
- the pull request will be reviewed

## Content

The manual is structured as followed:

```@contents
Pages = [
    "installation.md",
    "usage.md",
    "examples.md",
    "functions.md",
    "references.md"
    ]
Depth = 2
```


