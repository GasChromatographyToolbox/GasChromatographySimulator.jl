# GasChromatographySimulator.jl

[![CI](https://github.com/JanLeppert/GasChromatographySimulator.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/JanLeppert/GasChromatographySimulator.jl/actions/workflows/ci.yml)
[![codecov.io](http://codecov.io/github/JanLeppert/GasChromatographySimulator.jl/coverage.svg?branch=main)](http://codecov.io/github/JanLeppert/GasChromatographySimulator.jl?branch=main)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://JanLeppert.github.io/GasChromatographySimulator.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://JanLeppert.github.io/GasChromatographySimulator.jl/dev)

A package for the simulation of gas chromatography (GC) with additional velocity gradients produced by:
- thermal gradients
- film thickness gradient
- column diameter gradient

## Installation

To install the (unregistered) package type:

```julia
julia> ] add https://github.com/JanLeppert/GasChromatographySimulator.jl
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

## Citation

```
@misc{GCSimulator,
  title = {GasChromatographySimulator.jl},
  author = {Jan Leppert},
  howpublished = {\url{https://github.com/JanLeppert/GasChromatographySimulator.jl}}
}
```
