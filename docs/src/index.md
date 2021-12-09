# GasChromatographySimulator.jl Documentation

```@contents
```

## Functions 

### Structures
```@docs
GasChromatographySimulator.System
```

### Helper
```@docs
GasChromatographySimulator.load_solute_database
```

### Physical Model
```@docs
GasChromatographySimulator.pressure(x, t, T_itp, pin_itp, pout_itp, L, d, gas; ng=false)
```

### Solving
```@docs
GasChromatographySimulator.simulate(par)
```

### Results
```@docs
GasChromatographySimulator.peaklist(sol, par)

GasChromatographySimulator.peaklist(sol, peak, par)
```

## Index

```@index
```
