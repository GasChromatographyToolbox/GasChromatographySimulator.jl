# GasChromatographySimulator.jl

[![DOI](https://joss.theoj.org/papers/10.21105/joss.04565/status.svg)](https://doi.org/10.21105/joss.04565)
[![DOI](https://zenodo.org/badge/421858896.svg)](https://zenodo.org/badge/latestdoi/421858896)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://JanLeppert.github.io/GasChromatographySimulator.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://JanLeppert.github.io/GasChromatographySimulator.jl/dev)
[![CI](https://github.com/JanLeppert/GasChromatographySimulator.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/JanLeppert/GasChromatographySimulator.jl/actions/workflows/ci.yml)
[![codecov.io](http://codecov.io/github/JanLeppert/GasChromatographySimulator.jl/coverage.svg?branch=main)](http://codecov.io/github/JanLeppert/GasChromatographySimulator.jl?branch=main)

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

## Citation

```
@article{Leppert2022, 
  title = {GasChromatographySimulator.jl},
  author = {Jan Leppert}, 
  journal = {Journal of Open Source Software},
  year = {2022}, 
  volume = {7}, 
  number = {76}, 
  pages = {4565},
  publisher = {The Open Journal}, 
  doi = {10.21105/joss.04565}, 
  url = {https://doi.org/10.21105/joss.04565}, 
}
```
