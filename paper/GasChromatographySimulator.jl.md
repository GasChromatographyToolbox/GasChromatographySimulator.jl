---
title: 'GasChromatographySimulator.jl'
tags:
	- Julia
	- gas chromatography
	- separation
	- 
authors:
	- name: Jan Leppert
		orcid: ...
		affiliation: 1
affiliations:
	- name: Institute of Nutritional and Food Sciences, University of Bonn
		index: 1
date: 19 February 2022
bibliography: paper.bib
---

# Summary
`GasChromatographySimulator.jl` is Julia package [@Julia] to simulate the separation of a number of substances in a (one-dimensional) gas chromatographic (GC) system with programmed temperature $T(t)$ and pressure $p_{in}(t)$ resp. $p_{out}(t)$. The package also allows for spatial changes of the diameter of the GC column $d(x)$, of the film thickness of the stationary phase $d_f(x)$ and of the temperature $T(x)$.

Gas chromatography is a method to separate a mixture of substances by injecting the mixture into a gas stream (mobile phase), which passes along a tube (called column, length $L$ and diameter $d$) coated with a stationary phase (film thickness $d_f$). The substances of the mixture interact with the stationary phase by partition between mobile and stationary phase, resulting in different velocities of the substances. At the end of the column the separated substances are registered in a detector resulting in a chromatogram (retention time and signal/peak width).  

The modeling of GC separations is used for the prediction of retention times and widths of the signals (the chromatogram) and is of interest for method development, especially in multidimensional GC [@Hou:2018, @Jaramillo:2020, @Gaida:2021]. While the presented package only models the one dimensional GC separation, it can be used as the base for a multi-dimensional GC separation.

# Statement of need
`GasChromatographySimulator.jl` provides an interface to define a GC system consisting of a column (length, diameter, film thickness, type of stationary and mobile phase), program (temperature and pressure program, optional thermal spatial change) and substance parameters. By providing the thermodynamic parameters for the interaction of substances and stationary phase, which can be estimated by isothermal GC measurements, the separation of any mixture of substances can be simulated. 

The model is based on solving the ODEs for migration of the substances through the GC system $t(x)$ and for the development of the temporal variance of the substance distribution $\tau^2(x,t)$, using the Julia package `DifferentialEquations.jl` [@DifferentialEquations]. 
$$
\frac{dt}{dx} = r(x,t)
$$
and
$$
\frac{d\tau^2}{dx} = H(x,t)r(x,t) + 2 \tau^2(x,t) \frac{\partial r(x,t)}{\partial t}
$$
with $r$ the inverse substance velocity ($r=1/u$) and $H$ the local plate height. The basic equations building the model are presented in an earlier publication [@Leppert:2020].

A collection of `Pluto.jl` notebooks [@Pluto] are made available together with this package to provide a simple to use user interface to setup and simulate arbitrary GC systems.

# (Features and Functionality)
main features of the software (in the text above)

## Typical workflow
short example of setting up a GC simulation

The typical definition of of GC separation consists of four steps (four type structures, for which different constructor functions exist):
	- System `GasChromatographySimulator.Systems`
	- Program `GasChromatographySimulator.Program`
	- Substances `GasChromatographySimulator.Substance`
	- Options `GasChromatographySimulator.Options`
These parameter subsets are combined by `GasChromatographySimulator.Parameters`.

The simulation of the GC separation defined in `par::GasChromatographySimulator.Parameters` is executed by
`GasChromatographySimulator.simulate(par)``

Results
- solution of `DifferentialEquations.jl` (note of `t` <-> `x`)
- peaklist
- plot chromatogram
- plot local solution

# Comparison with existing software
Some simulations of GC exist [@Boswell:2012, @EZGC, @Gaida:2021, @Hou:2018] with different focus. 

An open-source software for the simulation of GC to predict retention times is the `GC Retention Predictor` [@Boswell:2012]. ... (JavaScript?, not easily to add other substances/stationary phases )

A propiatary software is the Pro EZGC Chromatogram Modeler from Restek [@ETZGC] (https://ez.restek.com/proezgc) ... (model is not known, a pre-defined set of stationary phases, from Restek, and substances ...)

Also in scientific papers specialized software is published (non-open source [@Gaida:2021], open source [@Hou:2018], both Matlab) ...

## (Conclusion)
like in every paper
# Acknowledgment
Jan Leppert is supported by the DFG research grant 452897652. (PRÜFEN)

## References
@Julia
@DifferentialEquations
@Pluto
@Boswell:2012 (weblink)
@EZGC
@Gaida:2021
@Hou:2018
@Leppert:2020

