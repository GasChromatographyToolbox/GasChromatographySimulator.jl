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

## Citation

@misc{GCSimulator,
  title = {GasChromatographySimulator.jl},
  author = {Jan Leppert},
  howpublished = {\url{https://github.com/JanLeppert/GasChromatographySimulator.jl}}
}