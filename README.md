# GasChromatographySimulator.jl

[![CI](https://github.com/JanLeppert/GasChromatographySimulator.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/JanLeppert/GasChromatographySimulator.jl/actions/workflows/ci.yml)
[![codecov.io](http://codecov.io/github/JanLeppert/GasChromatographySimulator.jl/coverage.svg?branch=main)](http://codecov.io/github/JanLeppert/GasChromatographySimulator.jl?branch=main)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://JanLeppert.github.io/GasChromatographySimulator.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://JanLeppert.github.io/GasChromatographySimulator.jl/dev)

A package for the simulation of gas chromatography (GC) with additional velocity gradients produced by:
- non-uniform temperature ``T(x)``
- non-uniform film thickness ``d_f(x)``
- non-uniform column diameter ``d(x)``

## Installation

To install the package type:

```julia
julia> ] add GasChromatographySimulator
```

To use the package type:

```julia
julia> using GasChromatographySimulator
```

## Documentation

Please read the [documentation page](https://janleppert.github.io/GasChromatographySimulator.jl/dev/) for more information.

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

## Citation

```
@misc{GCSimulator,
  title = {GasChromatographySimulator.jl},
  author = {Jan Leppert},
  howpublished = {\url{https://github.com/JanLeppert/GasChromatographySimulator.jl}}
}
```