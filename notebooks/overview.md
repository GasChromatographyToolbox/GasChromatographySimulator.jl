# Overview notebooks

- [x] `simulation_example.jl` - general example of simulation of a GC-system with optional thermal gradient (exponential/linear model)
- [x] `simulation_conventional_GC.jl` - simulation of a GC-system without thermal gradient and outlet pressure as "vacuum" or "atmospheric" -> [] base of a future notebook for teaching (there: definition of two GC-systems, one set of substances and compare the simulation results)
- [x] `simulation_conventional_GC_load_db.jl` - simulation of a GC-system without thermal gradient and loading of a self-defined substance database -> [] add the option of also loading measured retention times for a given GC-system
- [x] `simulation_conventional_GC_load_2dbs.jl` - simulation of a GC-system without thermal gradient and loading of two different substance databases and simulation of the common substances with the same GC-system and comparing the result -> [] add the option of also loading measured retention times for a given GC-system 
- [] `simulation_example_input_gradient_function.jl` - simulation of a GC-system with optional thermal gradient (self-defined equation)
- [] `test_abstol_reltol.jl` - test of the simulation with different absolute and relative tolerances for the solver of the ODEs -> [] copy this to another place?

For all notebooks:
- [] change the name to easier/shorter description