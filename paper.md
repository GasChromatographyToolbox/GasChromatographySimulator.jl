---
title: 'GasChromatographySimulator.jl'
tags:
  - Julia
  - gas chromatography
  - separation
  - thermodynamic parameters
  - simulation
authors:
  - name: Jan Leppert
    orcid: 0000-0001-8857-8103
    affiliation: 1
affiliations:
 - name: Institute of Nutritional and Food Sciences, University of Bonn
   index: 1
date: 28 March 2022
bibliography: paper.bib
---

# Summary
`GasChromatographySimulator.jl` is a Julia package [@Julia] to simulate the separation of a number of substances in a (one-dimensional) gas chromatographic (GC) system with programmed temperature $T(t)$ and programmed inlet pressure $p_{in}(t)$ (pressure controlled) resp. programmed flow $F(t)$ (flow controlled). In principle, the outlet pressure can also be programmed $p_{out}(t)$. The package also allows for spatial changes of the diameter of the GC column $d(x)$, of the film thickness of the stationary phase $d_f(x)$ and of the temperature $T(x)$. The simulation is based on solving a system of ordinary differential equations (ODE) for the migration of the substances through the GC system $t(x)$ and for the development of the temporal variance of the substance distribution $\tau^2(x,t)$

Gas chromatography is used in analytical chemistry as a method to separate a mixture of substances by injecting the mixture into a gas stream (mobile phase), which passes along a tube (called column, length $L$ and diameter $d$) coated with a stationary phase (film thickness $d_f$). The substances of the mixture interact with the stationary phase by partition between mobile and stationary phase, resulting in different velocities of the substances. At the end of the column the separated substances are registered in a detector resulting in a chromatogram (retention time and signal/peak width). Gas chromatography is used for qualitative and quantitative analysis, e.g. in petro chemistry, food chemistry, environmental analytics and forensic science.

The modeling of GC separations is used for the prediction of retention times and widths of the signals. With these results the separation of substances can be evaluated via the resolution of neighboring substance peaks in the chromatogram.  Such simulations are of interest for method development, especially in multidimensional GC [@Hou:2018; @Jaramillo:2020; @Gaida:2021]. While the presented package only models the one dimensional GC separation, it can be used as the base for a multi-dimensional GC separation.

# Statement of need
`GasChromatographySimulator.jl` provides an interface to define a GC system consisting of: 

- column (length, diameter, film thickness, type of stationary and mobile phase) 
- program (temperature and pressure program, optional thermal spatial change)
- substance parameters (thermodynamic parameters, diffusion coefficient)
- additional options (e.g. tolerances, algorithm for solving ODEs, model of viscosity) 

By providing the thermodynamic parameters [@Blumberg:2017] for the interaction of substances and stationary phase, which can be estimated by isothermal GC measurements, the separation of mixtures of substances on can be simulated for a wide range of GC systems. 

The ODE system for migration $t(x)$ and temporal variance development $\tau^2(x,t)$ is:  
$$
\frac{dt}{dx} = r(x,t)
$$
and
$$
\frac{d\tau^2}{dx} = H(x,t)r(x,t) + 2 \tau^2(x,t) \frac{\partial r(x,t)}{\partial t}
$$
with $r$ the inverse substance velocity ($r=1/u$) and $H$ the local plate height, [@Leppert:2020a], on the interval of $0 \leq x \leq L$, where $L$ is the length of the column. The basic equations building the model are presented in an earlier publication [@Leppert:2020b] and can be found in the documentation of the package. This ODE system is solved by using the Julia package `DifferentialEquations.jl` [@DifferentialEquations].

A collection of `Pluto.jl` notebooks [@Pluto] are made available together with this package to provide a simple to use user interface to setup and simulate arbitrary GC systems.

# Acknowledgment
Jan Leppert is supported by the DFG research grant 452897652.

# References






