# Overview notebooks

- `simulation_conventional_GC.jl` - simulation of a conventional GC-system (constant temperature, constant diameter and constant film thickness along the column) and outlet pressure as "vacuum" or "atmospheric"
- `simulation_conventional_GC_load_2dbs.jl` - simulation of a conventional GC-system (constant temperature, constant diameter and constant film thickness along the column) and loading of up to two different substance databases and simulation of the common substances with the same GC-system and comparing the result. Also, an option is given, to load measured retention times and compare these to the simulations.
- `simulation_example.jl` - general example of simulation of a GC-system with optional thermal gradient (exponential/linear model of temperature change along the column) and constant diameter and constant film thickness along the column. 
- `simulation_example_input_gradient_function.jl` - simulation of a GC-system with optional thermal gradient where the temperature change along the column is defined by a user-defined equation (cosine-function as example)
- 
- [] `test_abstol_reltol.jl` - test of the simulation with different absolute and relative tolerances for the solver of the ODEs -> [] copy this to another place?

