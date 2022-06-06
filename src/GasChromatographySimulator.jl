module GasChromatographySimulator

using Reexport
@reexport using Interpolations
using QuadGK
@reexport using DifferentialEquations
using ForwardDiff
@reexport using DataFrames
@reexport using CSV
using Plots
using HypertextLiteral
using PlutoUI
using ChemicalIdentifiers

# some constants
Tst = 273.15            # K
R = 8.31446261815324    # J mol⁻¹ K⁻¹
Tn = 25.0 + Tst         # K
pn = 101300             # Pa

# add a custom database for ChemicalIdentifiers 
#custom_database_url = "https://github.com/JanLeppert/GasChromatographySimulator.jl/blob/main/data/custom_CI_db.tsv"
filepath = string(pwd(), "/data/custom_CI_db.tsv")
#ChemicalIdentifiers.load_data!(:custom, url = custom_database_url)
ChemicalIdentifiers.load_data!(:custom, file = filepath)
ChemicalIdentifiers.load_db!(:custom)

# ---Begin-Structures---
"""
    Column(L, d, a_d, df, a_df, sp, gas)

Structure describing the GC Column. 

# Arguments

* `L`: Length of the capillary measured in m (meter)
* `d`: A function `d(x, a_d)` of `x`, the position along the capillary, describing the diameter in m (meter).
* `a_d`: Parameters of the diameter function. 
* `d_f`: A function `d_f(x, a_df)` of `x`, describing the film thickness in m (meter).
* `sp`: The name of the stationary phase.
* `gas`: The name of the mobile phase. Allowed values: He, H2 or N2.
"""
struct Column{Fd<:Function, Fdf<:Function}
    L                       # column length in m
    d::Fd                   # column internal diameter in m as function of x
    a_d::Array{Float64,1}   # parameters of the diameters function d(x)
    df::Fdf                 # column film thickness in m as function of x
    a_df::Array{Float64}    # parameters of the film thickness function df(x)
    sp::String              # stationary phase of the column
    gas::String             # gas of the mobile phase ["He", "H2", "N2"]
    Column(L, d, a_d, df, a_df, sp, gas) = (gas!="He" && gas!="H2" && gas!="N2") ? error("Wrong selection for 'gas'. Choose 'He', 'H2' or 'N2'.") : new{typeof(d), typeof(df)}(L, d, a_d, df, a_df, sp, gas)
end

"""
    Program(time_steps, temp_steps, Fpin_steps, pout_steps, gf, a_gf, T_itp, Fpin_itp, pout_itp)

Structure to describe the temperature and flow/pressure program of a GC Column. The function `gf` describes an optional thermal gradient.

# Arguments
* `time_steps`: Time steps in s (seconds). 
* `temp_steps`: Temperature steps in °C (degree celsius).
* `Fpin_steps`: Flow steps in m³/s resp. inlet pressure steps in Pa(a) (pascal, absolute).
* `pout_steps`: Outlet pressure steps in Pa(a) (pascal, absolute).
* `gf`: Gradient function `gf(x, a_gf)`, describes the thermal gradient.
* `a_gf`: Parameters of the gradient function.
* `T_itp`: Interpolated (linear) temperature `T(x,t)`, constructed from `time_steps`, `temp_steps` and `gf`.
* `Fpin_itp`: Interpolated (linear) flow/inlet pressure `Fpin(t)`, constructed from `time_steps` and `Fpin_steps`.
* `pout_itp`: Interpolated (linear) outlet pressure `pout(t)`, constructed from `time_steps` and `pout_steps`.  

Note: The length of the arrays `time_steps`, `temp_steps`, `pin_steps` and `a_gf`
have to be the same.
"""
struct Program{Fgf<:Function}
    time_steps::Array{<:Real, 1}            # vector time steps for the temperature program
    temp_steps::Array{<:Real, 1}            # vector temperature steps for the temperature program 
    Fpin_steps::Array{<:Real, 1}            # vector flow resp. inlet pressure steps for the pressure program
    pout_steps::Array{<:Real, 1}            # vector outlet pressure steps for the pressure program
    gf::Fgf                                 # function of x of the gradient form
    a_gf::Array{<:Real}                     # parameters of the gradient function gf(x)
  	T_itp::Interpolations.Extrapolation     # interpolation function of T(x,t)
    Fpin_itp::Interpolations.Extrapolation  # interpolation function of Fpin(t)
    pout_itp::Interpolations.Extrapolation	# interpolation function of pout(t)
    # add inner constructor to check the lengths of the Arrays and of the result of gf
    Program(ts, Ts, Fpis, pos, gf, agf, T, Fpin, po) = (length(ts)!=length(Ts) || length(ts)!=length(gf(0.0)) || length(ts)!=length(Fpis) || length(ts)!=length(pos)) || length(ts)!=size(agf)[1] ? error("Mismatch between length(time_steps) = $(length(ts)), length(temp_steps) = $(length(Ts)), length(Fpin_steps) = $(length(Fpis)), length(pout_steps) = $(length(pos)), length(gf(0.0)) = $(length(gf(0.0))) and size(a_gf)[1] = $(size(agf)[1]).") : new{typeof(gf)}(ts, Ts, Fpis, pos, gf, agf, T, Fpin, po)
end

"""
    Substance(name, CAS, Tchar, θchar, ΔCp, φ₀, ann, Dag, t₀, τ₀)

Structure to describe the properties of a solute, which migrates through the GC Column. These datas are in most cases read from a database with the function `load_solute_database()`.

# Arguments
* `name`: Name of the solute. 
* `CAS`: CAS number of the solute.
* `Tchar`: Characterisic temperature (in K). One of the three distribution-centric thermodynamic parameters describing the retention of this solute on the given stationary phase.
* `θchar`: Characterisic parameters (in °C). One of the three distribution-centric thermodynamic parameters describing the retention of this solute on the given stationary phase.
* `ΔCp`: Change of the isobaric heat capacity moving from the mobile to the stationary phase (in J mol⁻¹ K⁻¹). One of the three distribution-centric thermodynamic parameters describing the retention of this solute on the given stationary phase.
* `φ₀`: Dimensionless film thickness (φ ≈ df/d) of the column for which the thermodynamic parameters (Tchar, θchar, ΔCp) were estimated.
* `ann`: Annotations. In most cases the source of the data is noted here.
* `Dag`: The diffusitivity of the solute `a` in the mobile phase `g` (in...). It is calculated by the function `diffusitivity()`.
* `t₀`: Initial time of the solute (in s) at the start of the simulation.
* `τ₀`: Initial peak width of the solute (in s) at the start of the simulation. 

See also: [`load_solute_database`](@ref)
"""
struct Substance
    name::String        # name of solute
    CAS::String         # CAS registry number
    Tchar               # characteristic temperature in K
    θchar               # characteristic thermal constant in °C
    ΔCp                 # 3rd parameter
    φ₀                  # dimless film thickness for which Tchar, θchar and ΔCp were estimated
    ann::String         # annotations, e.g. the source of the data from which Tchar, θchar and ΔCp were estimated
    Dag                 # diffusion coefficient of analyt in a gas, calculate from structure (or from measurements)
    t₀                  # initial time in s  	
    τ₀                  # initial peak width in s   
end

"""
    Options(alg, abstol, reltol, Tcontrol, odesys, ng, vis)

Structure describing some general options for the simulation. 

# Arguments
* `alg`: The algorithm used for the ODE solver. The algorithms `OwrenZen3()`, `OwrenZen4()` and `OwrenZen5()` are recommended.
* `abstol`: The absolute tolerance for the ODE solver. Recommended value 1e-6 to 1e-8.
* `reltol`: The relative tolerance for the ODE solver. Recommended value 1e-3 to 1e-5. 
* `Tcontrol`: Option defining at which point of the column the temperature program is calculated. The options are `inlet` (x=0) and `outlet` (x=L).
* `odesys`: Combine the ODEs for migration and peak-width into a system of ODEs (`odesys = true`) or solve the two ODEs separately (`odesys = false`).
* `ng`: Option to calculate the simulation without a gradient (`ng = true`) or with a gradient (`ng = false`). This distinction is made because of partly manuall differentiation (problem of automatic differentiation with integrals, e.g. in the `flow_restriction()` function. -> **TODO**: test package Quadrature.jl as alternative to QuadGK.jl for integration)
* `vis`: Used model of viscosity. `HP` is a second-order polynomial taken from the HP flow calculator. `Blumberg` is an emperical formula according to the book
    `Temperature-programmed Gas Chromatography` by Leonid M. Blumberg (2010, Wiley-VCH) 
* `control`: Control of the "Flow" or of the "Pressure" (at column inlet) during the program

**TODO**: add option for the retention model ('ABC', 'K-centric')

For more informations about the arguments `alg`, `abstol` and `reltol` see the documentation of the DifferentialEquations.jl package.
"""
struct Options
    alg                 # algorithmen for the ODE solver
    abstol              # absolute tolerance for ODE solver
    reltol              # relative tolerance for ODE solver 
    Tcontrol::String    # temperature control at 'inlet' (top) or 'outlet' (bottom) of the column
	odesys::Bool  		# calculate the two ODEs (migration and peak-width) separately (false) or 
                        # combined as a system of ODEs (true)                        
    ng::Bool            # non-gradient calculation, ignores a defined spatial change of d, df or T
    vis::String         # viscosity model 'HP' or 'Blumberg'
    control::String     # control of the 'Flow' or of the inlet 'Pressure' during the program
end

"""
    Parameters(col, prog, sub, opt)

Structure describing all parameters for the simulation of a GC system. 

# Arguments
* `col`: Structure `Column` describing the parameters of the GC column and
    mobile phase gas.
* `prog`: Structure `Program` describing the temperature and pressure
    program of a GC Column.
* `sub`: An array of the structure `Substance` describing the parameters of
    the solutes which are separated in the GC simulation. 
* `opt`: Structure `Options` describing additional option parameters.
"""
struct Parameters
    col::Column
    prog::Program
    sub::Array{Substance,1}
    opt::Options
end
# ---End-Structures---

# ---Begin-Constructor-functions-for-Structures---

"""
    Column(L, d, df, sp, gas)

Construct the structure `Column` with given values for the case
of constant diameter `d` and film thickness `df`. 

# Arguments
* `L`: Length of the capillary measured in m (meter).
* `d`: Diameter of the capillary measured in m (meter). 
* `d_f`: Film thickness of the capillary measured in m (meter).
* `sp`: The name of the stationary phase.
* `gas`: The name of the mobile phase. Allowed values: He, H2 or N2.

# Examples
```julia
julia> Column(10.0, 0.1e-3, 0.1e-6, "DB5", "He")
	```
"""
function Column(L, d, df, sp, gas)
    # function to construct the Column structure
    # for the case of constant diameter and constant film thickness
    d_func(x) = gradient(x, [d])
    df_func(x) = gradient(x, [df])
    col = Column(L, d_func, [d], df_func, [df], sp, gas)
    return col
end

"""
    Program(time_steps, temp_steps, Fpin_steps, pout_steps, a_gf, Tcontrol, L)

Construct the structure `Program` with given values. 

# Arguments
* `time_steps`: Time steps in s (seconds). 
* `temp_steps`: Temperature steps in °C (degree celsius).
* `Fpin_steps`: Flow steps in m³/s resp. inlet pressure steps in Pa(a).
* `pout_steps`: Outlet pressure steps in Pa(a) (pascal, absolute).
* `a_gf`: Parameters of the gradient function.
* `Tcontrol`: Option defining at which point of the column the temperature
program is calculated. The options are `inlet` (x=0) and `outlet` (x=L).
* `L`: Length of the capillary measured in m (meter).

The length of the arrays `time_steps`, `temp_steps`, `pin_steps`, `pout_steps` and `a_gf`
have to be the same.

The arguments `Tcontrol` and `L` are used to construct the thermal gradient
function `gf(x)` and the temperature interpolation `T_itp(x,t)`.

# Examples
```julia
julia> Program([0.0, 60.0, 300.0, 120.0],
        [40.0, 40.0, 320.0, 320.0],
        300000.0.*ones(4),
        zeros(4),
        [[20.0, 20.0, 20.0, 20.0] zeros(4) 10.0.*ones(4) [0.0, -2.0, -5.0, -5.0]],
        "inlet",
        10.0)
```
"""
function Program(time_steps::Array{<:Real, 1}, temp_steps::Array{<:Real, 1}, Fpin_steps::Array{<:Real, 1}, pout_steps::Array{<:Real, 1}, a_gf::Array{<:Real, 2}, Tcontrol, L)
    # function to construct the Program structure
    # using as gradient function the exponential model 'gf_exp(x,a_gf,Tcontrol)'
    gf(x) = gradient(x, a_gf; Tcontrol=Tcontrol)
    T_itp = temperature_interpolation(time_steps, temp_steps, gf, L)
    Fpin_itp = steps_interpolation(time_steps, Fpin_steps)
    pout_itp = steps_interpolation(time_steps, pout_steps)
    prog = Program(time_steps, temp_steps, Fpin_steps, pout_steps, gf, a_gf, T_itp, Fpin_itp, pout_itp)
    return prog
end

"""
    Program(time_steps, temp_steps, Fpin_steps, pout_steps, ΔT_steps, x₀_steps, L₀_steps, α_steps, Tcontrol, L) 

Construct the structure `Program` with given values. 

# Arguments
* `time_steps`: Time steps in s (seconds). 
* `temp_steps`: Temperature steps in °C (degree celsius).
* `Fpin_steps`: Flow steps in m³/s resp. inlet pressure steps in Pa(a).
* `pout_steps`: Outlet pressure steps in Pa(a) (pascal, absolute).
* `ΔT_steps`: Parameters of the gradient function. Temperature difference in °C.
* `x₀_steps`: Parameters of the gradient function. Spatial offset of the gradient in m.
* `L₀_steps`: Parameters of the gradient function. Distance over which the temperature difference in ΔT_steps is measured, in m.
* `α_steps`: Parameters of the gradient function. Factor describing the gradient profile. 
* `Tcontrol`: Option defining at which point of the column the temperature
program is calculated. The options are `inlet` (x=0) and `outlet` (x=L).
* `L`: Length of the capillary measured in m (meter).

The length of the arrays `time_steps`, `temp_steps`, `Fpin_steps`, `pout_steps`, `ΔT_steps`, `x₀_steps`, `L₀_steps` and `α_steps`
have to be the same.

The arguments `Tcontrol` and `L` are used to construct the thermal gradient
function `gf(x)` and the temperature interpolation `T_itp(x,t)`.

# Examples
```julia
julia> Program([0.0, 60.0, 300.0, 120.0],
        [40.0, 40.0, 320.0, 320.0],
        300000.0.*ones(4),
        zeros(4),
        [20.0, 20.0, 20.0, 20.0],
        zeros(4),
        10.0.*ones(4),
        [0.0, -2.0, -5.0, -5.0],
        "inlet",
        10.0)
```
"""
function Program(time_steps::Array{<:Real, 1}, temp_steps::Array{<:Real, 1}, Fpin_steps::Array{<:Real, 1}, pout_steps::Array{<:Real, 1}, ΔT_steps::Array{<:Real, 1}, x₀_steps::Array{<:Real, 1}, L₀_steps::Array{<:Real, 1}, α_steps::Array{<:Real, 1}, Tcontrol, L)
    # function to construct the Program structure
    # using as gradient function the exponential model 'gf_exp(x,a_gf,Tcontrol)'
    a_gf = [ΔT_steps x₀_steps L₀_steps α_steps]
    gf(x) = gradient(x, a_gf; Tcontrol=Tcontrol)
    T_itp = temperature_interpolation(time_steps, temp_steps, gf, L)
    Fpin_itp = steps_interpolation(time_steps, Fpin_steps)
    pout_itp = steps_interpolation(time_steps, pout_steps)
    prog = Program(time_steps, temp_steps, Fpin_steps, pout_steps, gf, a_gf, T_itp, Fpin_itp, pout_itp)
    return prog
end

"""
    Program(TP, FpinP, poutP, ΔTP, x₀P, L₀P, αP, Tcontrol, L; time_unit="min") 

Construct the structure `Program` with conventional formulation (see [`conventional_program`](@ref)) of programs. 

# Arguments
* `TP`: conventional formulation of a temperature program. 
* `FpinP`: conventional formulation of a Flow (in m³/s) resp. inlet pressure (in Pa(a)) program.
* `poutP`: conventional formulation of a outlet pressure (in Pa(a)) program.
* `ΔTP`: conventional formulation of a program for temperature gradient parameter ΔT. 
* `x₀P`: conventional formulation of a program for temperature gradient parameter x₀.
* `L₀P`: conventional formulation of a program for temperature gradient parameter L₀.
* `αP`: conventional formulation of a program for temperature gradient parameter αP.
* `Tcontrol`: Option defining at which point of the column the temperature
program is calculated. The options are `inlet` (x=0) and `outlet` (x=L).
* `L`: Length of the capillary measured in m (meter).
* `time_unit`: unit of time in the programs, "min"` (default) times are measured in minutes, "s" times are measured in seconds.


The arguments `Tcontrol` and `L` are used to construct the thermal gradient
function `gf(x)` and the temperature interpolation `T_itp(x,t)`.

"""
function Program(TP::Array{<:Real, 1}, FpinP::Array{<:Real, 1}, poutP::Array{<:Real, 1}, ΔTP::Array{<:Real, 1}, x₀P::Array{<:Real, 1}, L₀P::Array{<:Real, 1}, αP::Array{<:Real, 1}, Tcontrol::String, L::Float64; time_unit="min")
    # using as gradient function the exponential model 'gf_exp(x,a_gf,Tcontrol)
    ts = Array{Array{Float64,1}}(undef, 7)
    ts[1], Ts = GasChromatographySimulator.conventional_program(TP; time_unit=time_unit)
    ts[2], Fps = GasChromatographySimulator.conventional_program(FpinP; time_unit=time_unit)
    ts[3], pouts = GasChromatographySimulator.conventional_program(poutP; time_unit=time_unit)
    ts[4], ΔTs = GasChromatographySimulator.conventional_program(ΔTP; time_unit=time_unit)
    ts[5], x₀s = GasChromatographySimulator.conventional_program(x₀P; time_unit=time_unit)
    ts[6], L₀s = GasChromatographySimulator.conventional_program(L₀P; time_unit=time_unit)
    ts[7], αs = GasChromatographySimulator.conventional_program(αP; time_unit=time_unit)
    time_steps = Float64[]
    for i=1:7
        time_steps = GasChromatographySimulator.common_time_steps(time_steps, ts[i])
    end
    temp_steps = GasChromatographySimulator.new_value_steps(Ts, ts[1], time_steps)
    Fpin_steps = GasChromatographySimulator.new_value_steps(Fps, ts[2], time_steps)
    pout_steps = GasChromatographySimulator.new_value_steps(pouts, ts[3], time_steps)
    ΔT_steps = GasChromatographySimulator.new_value_steps(ΔTs, ts[4], time_steps)
    x₀_steps = GasChromatographySimulator.new_value_steps(x₀s, ts[5], time_steps)
    L₀_steps = GasChromatographySimulator.new_value_steps(L₀s, ts[6], time_steps)
    α_steps = GasChromatographySimulator.new_value_steps(αs, ts[7], time_steps)
   
    a_gf = [ΔT_steps x₀_steps L₀_steps α_steps]
    gf(x) = gradient(x, a_gf; Tcontrol=Tcontrol)
    T_itp = temperature_interpolation(time_steps, temp_steps, gf, L)
    Fpin_itp = steps_interpolation(time_steps, Fpin_steps)
    pout_itp = steps_interpolation(time_steps, pout_steps)
    prog = Program(time_steps, temp_steps, Fpin_steps, pout_steps, gf, a_gf, T_itp, Fpin_itp, pout_itp)
    return prog
end

"""
    Program(time_steps::Array{<:Real, 1}, temp_steps::Array{<:Real, 1}, Fpin_steps::Array{<:Real, 1}, pout_steps::Array{<:Real, 1}, L)

Construct the structure `Program` with given values for the case
without a thermal gradient. 

# Arguments
* `time_steps::Array{<:Real, 1}`: Time steps in s (seconds). 
* `temp_steps::Array{<:Real, 1}`: Temperature steps in °C (degree celsius).
* `Fpin_steps::Array{<:Real, 1}`: Flow steps in m³/s resp. inlet pressure steps in Pa(a).
* `pout_steps::Array{<:Real, 1}`: Outlet pressure steps in Pa(a) (pascal, absolute).
* `L`: Length of the capillary measured in m (meter).

The length of the arrays `time_steps`, `temp_steps`, `Fpin_steps` and `pout_steps`
have to be the same.

The argument `L` is used to construct the temperature interpolation `T_itp(x,t)`.

# Examples
```julia
julia> Program([0.0, 60.0, 300.0, 120.0],
        [40.0, 40.0, 320.0, 320.0],
        300000.0.*ones(4),
        zeros(4),
        10.0)
```
"""
function Program(time_steps::Array{<:Real, 1}, temp_steps::Array{<:Real, 1}, Fpin_steps::Array{<:Real, 1}, pout_steps::Array{<:Real, 1}, L)
    # function to construct the Program structure
    # without a thermal gradient
    # using as gradient function the exponential model 'gf_exp(x,a_gf,Tcontrol)'
    a_gf = [zeros(length(time_steps)) zeros(length(time_steps)) L.*ones(length(time_steps)) zeros(length(time_steps))]
    gf(x) = gradient(x, a_gf)
    T_itp = temperature_interpolation(time_steps, temp_steps, gf, L)
    Fpin_itp = steps_interpolation(time_steps, Fpin_steps)
    pout_itp = steps_interpolation(time_steps, pout_steps)
    prog = Program(time_steps, temp_steps, Fpin_steps, pout_steps, gf, a_gf, T_itp, Fpin_itp, pout_itp)
    return prog
end

"""
    Program(TP, FpinP, L; pout="vacuum", time_unit="min")

Construct the structure `Program` with conventional formulation (see [`conventional_program`](@ref)) of programs for the case
without a thermal gradient. 

# Arguments
* `TP`: conventional formulation of a temperature program. 
* `FpinP`: conventional formulation of a Flow (in m³/s) resp. inlet pressure (in Pa(a)) program.
* `L`: Length of the capillary measured in m (meter).
* `pout`: Outlet pressure, "vacuum" (default) or "atmosphere".
* `time_unit`: unit of time in the programs, "min"` (default) times are measured in minutes, "s" times are measured in seconds.

The argument `L` is used to construct the temperature interpolation `T_itp(x,t)`.

# Examples
```julia
julia> Program((40.0, 1.0, 5.0, 280.0, 2.0, 20.0, 320.0, 2.0),
                (400000.0, 10.0, 5000.0, 500000.0, 20.0),
                10.0)
```
"""
function Program(TP, FpinP, L; pout="vacuum", time_unit="min")
    ts1, Ts = GasChromatographySimulator.conventional_program(TP; time_unit=time_unit)
    ts2, Fps = GasChromatographySimulator.conventional_program(FpinP; time_unit=time_unit)
    time_steps = GasChromatographySimulator.common_time_steps(ts1, ts2)
    temp_steps = GasChromatographySimulator.new_value_steps(Ts, ts1, time_steps)
    Fpin_steps = GasChromatographySimulator.new_value_steps(Fps, ts2, time_steps)
    if pout == "vacuum"
        pout_steps = zeros(length(time_steps))
    else
        pout_steps = pn.*ones(length(time_steps))
    end
    a_gf = [zeros(length(time_steps)) zeros(length(time_steps)) L.*ones(length(time_steps)) zeros(length(time_steps))]
    gf(x) = gradient(x, a_gf)
    T_itp = temperature_interpolation(time_steps, temp_steps, gf, L)
    Fpin_itp = steps_interpolation(time_steps, Fpin_steps)
    pout_itp = steps_interpolation(time_steps, pout_steps)
    prog = Program(time_steps, temp_steps, Fpin_steps, pout_steps, gf, a_gf, T_itp, Fpin_itp, pout_itp)
    return prog
end

"""
    Options(;alg=OwrenZen5(), abstol=1e-6, reltol=1e-3, Tcontrol="inlet", odesys=true, ng=false, vis="Blumberg", control="Pressure")

Construct the structure `Options` with default values. 

# Arguments
* `alg`: The algorithm used for the ODE solver. The algorithms
    `OwrenZen3()`, `OwrenZen4()` and `OwrenZen5()` are recommended.
* `abstol`: The absolute tolerance for the ODE solver. Recommended value
    1e-6 to 1e-8.
* `reltol`: The relative tolerance for the ODE solver. Recommended value
1e-3 to 1e-5. 
* `Tcontrol`: Option defining at which point of the column the temperature
    program is calculated. The options are `inlet` (x=0) and `outlet` (x=L).
* `odesys`: Combine the ODEs for migration and peak-width into a system of
    ODEs (`odesys = true`) or solve the two ODEs separately (`odesys = false`).
* `ng`: Option to calculate the simulation without a gradient (`ng = true`)
    or with a gradient (`ng = false`). This distinction is made because of
    partly manuall differentiation (problem of automatic differentiation with
    integrals, e.g. in the `flow_restriction()` function. -> TODO: test
    package Quadrature.jl as alternative to QuadGK.jl for integration)
* `vis`: Used model of viscosity. `HP` is a second-order polynomial taken from the HP flow calculator. `Blumberg` is an emperical formula according to the book
    `Temperature-programmed Gas Chromatography` by Leonid M. Blumberg (2010, Wiley-VCH) 
* `control`: Control of the "Flow" or of the "Pressure" (at column inlet) during the program

For more informations about the arguments `alg`, `abstol` and `reltol` see
the documentation of the DifferentialEquations.jl package.

# Examples
```julia
julia> Options()
```

```julia
julia> Options(abstol=1e-8, Tcontrol="outlet")
```
"""
function Options(;alg=OwrenZen5(), abstol=1e-6, reltol=1e-3, Tcontrol="inlet", odesys=true, ng=false, vis="Blumberg", control="Pressure")
    opt = Options(alg, abstol, reltol, Tcontrol, odesys, ng, vis, control)
    return opt
end

"""
    Options(alg, abstol, reltol, Tcontrol, odesys; ng=false, vis="Blumberg", control="Pressure")

Construct the structure `Options` with given values. 

# Arguments
* `alg`: The algorithm used for the ODE solver. The algorithms
    `OwrenZen3()`, `OwrenZen4()` and `OwrenZen5()` are recommended.
* `abstol`: The absolute tolerance for the ODE solver. Recommended value
    1e-6 to 1e-8.
* `reltol`: The relative tolerance for the ODE solver. Recommended value
1e-3 to 1e-5. 
* `Tcontrol`: Option defining at which point of the column the temperature
    program is calculated. The options are `inlet` (x=0) and `outlet` (x=L).
* `odesys`: Combine the ODEs for migration and peak-width into a system of
    ODEs (`odesys = true`) or solve the two ODEs separately (`odesys = false`).
* `ng`: Option to calculate the simulation without a gradient (`ng = true`)
    or with a gradient (`ng = false`). This distinction is made because of
    partly manuall differentiation (problem of automatic differentiation with
    integrals, e.g. in the `flow_restriction()` function. -> TODO: test
    package Quadrature.jl as alternative to QuadGK.jl for integration)
* `vis`: Used model of viscosity. `HP` is a second-order polynomial taken from the HP flow calculator. `Blumberg` is an emperical formula according to the book
    `Temperature-programmed Gas Chromatography` by Leonid M. Blumberg (2010, Wiley-VCH) 
* `control`: Control of the "Flow" or of the "Pressure" (at column inlet) during the program

For more informations about the arguments `alg`, `abstol` and `reltol` see
the documentation of the DifferentialEquations.jl package.

# Examples
```julia
julia> Options(OwrenZen3(), 1e-7, 1e-4, "inlet", true)
```

```julia
julia> Options(OwrenZen3(), 1e-7, 1e-4, "inlet", true; ng=true, vis="HP", control="Flow")
```
"""
function Options(alg, abstol, reltol, Tcontrol, odesys; ng=false, vis="Blumberg", control="Pressure")
    opt = Options(alg, abstol, reltol, Tcontrol, odesys, ng, vis, control)
    return opt
end

# Aliases for compatibility
constructor_System = Column
constructor_Program = Program
constructor_Options = Options
#---End-Constructor-functions-for-Structures---

#---Begin-Functions-used-for-Parameter-construction 
"""
    gradient(x, a; Tcontrol="inlet")

Defines the gradient of the column diameter, film thickness or
temperature along the GC column.  

# Arguments
* `x`: Position along the GC column, in m.
* `a`: Parameters of the gradient function

The form of `a` decides the actual used function for the gradient. The
following options are realized:
* `a` is a single value (e.g. Float or Int): The gradient function is constant for all positions `x` with the value of `a`. 
* `a` is a 1D-array of length = 4: The gradient function is a exponential function. The 4 values of `a`: `a = [f₀, x₀, L₀, α]` with 
    * `f₀`: Start value at `x = x₀`.
    * `x₀`: Shift in `x` position.
    * `L₀`: Distance over which the value drops from f₀ to 0.
    * `α`: Factor describing the gradient profile.
    if `α <= 0` 
        `f = f₀ .* (1 .- exp.(α.*(1 .- (x.-x₀) ./ L₀)) + (1 .- (x.-x₀) ./ L₀)
        .* exp.(α))`
    if `α > 0`
        `f = f₀ .* (exp.(-α.*(x.-x₀) ./ L₀) .- (x.-x₀) ./ L₀ .* exp.(-α))` 
* `a` is a 2D-array with size = (n, 4): The gradient function is an
    exponential function. The 4 entrys have the same meaning as above, but their
    values can change over the times defined in the `time_steps` of the
    Program structure. At these `time_steps[i]` the gradient function is
    described by the corresponding parameters `a[i,:]`. Between the
    `time_steps[i]` and `time_steps[j]` the value of the gradient function at position `x` is linear
    interpolated from the gradien functions defined by `a[i,:]` and `a[j,:]`.

# Examples
```julia
julia> d(x) = gradient(x, 0.1e-3)
```

```julia
julia> gf(x) = gradient(x, [[20.0, 20.0, 20.0, 20.0] zeros(4) 10.0.*ones(4) [0.0, -2.0, -5.0, -5.0]])
```
"""
# gradient functions
function gradient(x, a; Tcontrol="inlet")
    if size(a) == (1,)
        # constant value, no gradient
        f = a[1]
    elseif size(a) > (1,)
        if  length(size(a))==1
            if length(a)==4
                # for diameter or film thickness, values of parameters 'a' are fixed
                # over time
                f₀ = a[1] # start value
                x₀ = a[2] # shift in x 
                L₀ = a[3] # length of the linear segment
                α  = a[4] # exponential factor, α=0 -> linear
                if α<=0 # concave/linear, weak change at beginning and strong change at end of column
                    f = f₀ .* (1 .- exp.(α.*(1 .- (x.-x₀) ./ L₀)) + (1 .- (x.-x₀) ./ L₀) .* exp.(α))
                elseif α>0 # convex function, strong change at beginning and weak change at end of the column
                    f = f₀ .* (exp.(-α.*(x.-x₀) ./ L₀) .- (x.-x₀) ./ L₀ .* exp.(-α))
                end
            # other functions ...
            end
        elseif length(size(a)) == 2
            if a==zeros(size(a))
                f = zeros(size(a)[1])
                #if size(a)[2] == 1
                # for thermal gradient, no change in time  of values of the
                # parameters
            elseif size(a)[2] == 4
                # for thermal gradient, values of parameters 'a[i,:]' can change
                # over time
                f₀ = a[:,1] # start value
                x₀ = a[:,2] # shift in x (e.g. to correct for the real position of IR-sensors)
                L₀ = a[:,3] # length of the linear segment
                α  = a[:,4] # exponential factor, α=0 -> linear
                f = Array{Float64}(undef, length(α))
                for i=1:length(α)
                    if α[i]<=0 # concave function, weak change at beginning and strong change at end of column
                        if Tcontrol=="inlet"
                            f[i] = f₀[i] .* (.- exp.(α[i].*(1 .- (x.-x₀[i]) ./ L₀[i])) + (1 .- (x.-x₀[i]) ./ L₀[i]) .* exp.(α[i]))
                        elseif Tcontrol=="outlet"
                            f[i] = f₀[i] .* (1 .- exp.(α[i].*(1 .- (x.-x₀[i]) ./ L₀[i])) + (1 .- (x.-x₀[i]) ./ L₀[i]) .* exp.(α[i]))
                        end
                    elseif α[i]>0 # convex function, strong change at beginning and weak change at end of the column
                        if Tcontrol=="inlet"
                            f[i] = f₀[i] .* (exp.(-α[i].*(x.-x₀[i]) ./ L₀[i]) .- (x.-x₀[i]) ./ L₀[i] .* exp.(-α[i]) .- 1)
                        elseif Tcontrol=="outlet"
                            f[i] = f₀[i] .* (exp.(-α[i].*(x.-x₀[i]) ./ L₀[i]) .- (x.-x₀[i]) ./ L₀[i] .* exp.(-α[i]))
                        end
                    end
                end
            # other functions ...
            end
        end
    end
    return f
end

"""
    temperature_interpolation(time_steps, temp_steps, gradient_function, L)

Construct the temperature function depending on position `x` and
time `t`.  

# Arguments
* `time_steps::Array{<:Real,1}`: Time steps in s (seconds). 
* `temp_steps::Array{<:Real,1}`: Temperature steps in °C (degree celsius).
* `gf::Function`: Gradient function `gf(x, a_gf)`, describes the thermal gradient.
* `L::Float64`: Length of the capillary measured in m (meter).

For the spatial dependency of the interpolated temperature `T_ipt(x,t)` the
gradient function `gf` is calculated every 1e-3 m (1 mm). Positions
inbetween are linear interpolated. For the temporal dependency the
temperatures `temp_steps` defined at the `time_steps` are linear
interpolated over time `t`.   

# Examples
```julia
julia> T_itp = temperature_interpolation([0.0, 60.0, 300.0, 120.0], [40.0, 40.0, 320.0, 320.0], gf, 10.0)
```
"""
function temperature_interpolation(time_steps::Array{<:Real,1}, temp_steps::Array{<:Real,1}, gradient_function::Function, L)
	T(x) = temp_steps .+ gradient_function(x) 
	nx = 0.0:1e-3:L # mm exact
	nt = cumsum(time_steps)
	Tmat = Array{Float64}(undef, length(nx), length(nt))
	for j=1:length(nt)
		for i=1:length(nx)
			Tmat[i,j] = T(nx[i])[j] + 273.15
		end
	end
	T_itp = LinearInterpolation((nx, nt), Tmat, extrapolation_bc=Flat())
	return T_itp
end

"""
    step_interpolation(time_steps, steps)

Construct a linear interpolated function depending on time `t` of the `steps`-values over `time_steps`.  

# Arguments
* `time_steps::Array{<:Real,1}`: Time steps in s (seconds). 
* `steps::Array{<:Real,1}`: steps, e.g. pressure or flow. 

# Examples
```julia
julia> pin_itp = steps_interpolation([0.0, 60.0, 300.0, 120.0], 
                                    [300000.0, 300000.0, 400000.0, 400000.0])
```
"""
function steps_interpolation(time_steps::Array{<:Real,1}, steps::Array{<:Real,1})
    s_itp = LinearInterpolation((cumsum(time_steps), ), steps, extrapolation_bc=Flat())
    return s_itp
end

"""
	CAS_identification(data::DataFrame)

Look up the substance name from the `data` dataframe with ChemicalIdentifiers.jl to find the `CAS`-number, the `formula`, the molecular weight `MW` and the `smiles`-identifier. If the name is not found in the database of ChemicalIdentifiers.jl a list with alternative names (`shortnames.csv`) is used. If there are still no matches, `missing` is used.
"""
function CAS_identification(Name::Array{String})
	shortnames = DataFrame(CSV.File("/Users/janleppert/Documents/GitHub/GasChromatographySimulator/data/shortnames.csv"))
	CAS = Array{Union{Missing,AbstractString}}(missing, length(Name))
	for i=1:length(Name)
		if Name[i] in shortnames.shortname
			j = findfirst(Name[i].==shortnames.shortname)
			ci = search_chemical(string(shortnames.name[j]))
		else
			ci = search_chemical(Name[i])
		end
		CAS[i] = string(ci.CAS[1], "-", ci.CAS[2], "-", ci.CAS[3])
	end
	id = DataFrame(Name=Name, CAS=CAS)
	return id
end

"""
    load_solute_database(db, sp, gas, solutes, t₀, τ₀)

Load the data of `solutes` for the stationary phase `sp` and the mobile
phase `gas` from the database `db` into an array of the structure `Substance`.

# Arguments
* `db::DataFrame`: DataFrame of the database. 
* `sp::String`: Name of the stationary phase.
* `gas::String`: Name of the mobile phase.
* `solutes::Array{<:AbstractString,1}`: Name of the solutes.
* `t₀::Array{Float64,1}`: Initial start times of the solutes.
* `τ₀::Array{Float64,1}`: Initial peak widths of the solutes. 

# Examples
```julia
julia> sub = load_solute_database(db, "DB5", "He", ["C10", "C11"], [0.0, 0.0], [0.5, 0.5])
```
"""
function load_solute_database(db::DataFrame, sp::String, gas::String, solutes::Array{<:AbstractString,1}, t₀::Array{Float64,1}, τ₀::Array{Float64,1})
	# compare names to CAS, using ChemicalIdentifiers.jl and a synonym list (different names of the same solute for different stationary phases) 
    id = CAS_identification(solutes)
    if sp == "" # no stationary phase is selected, like for transferlines
        Name = id.Name
        CAS = id.CAS
        # use placeholder values
        Tchar = ones(length(Name))
        θchar = ones(length(Name))
        ΔCp = ones(length(Name))
        φ₀ = ones(length(Name))
        Annotation = fill("no sp", length(Name))
        t₀_ = t₀
        τ₀_ = τ₀ 
    elseif size(db)[2]>14
        error("Data format not supported. Use the appended database structure.")
	elseif isa(findfirst(unique(db.Phase).==sp), Nothing) && sp!=""
		error("Unknown selction of stationary phase. Choose one of these: $phases_db")	
	else
        # old database config (simpler, additional stationary phases are appended in new rows, multiple entrys for a solute)
        # 1. Filter the stationary phase
        db_filtered = filter([:Phase] => x -> x==sp, db)
        # 2. Filter the solutes
        db_filtered_1=db_filtered[in(id.CAS).(db_filtered.CAS),:]
        # values
        Name = db_filtered_1.Name
        CAS = db_filtered_1.CAS
        Tchar = db_filtered_1.Tchar.+Tst
        θchar = db_filtered_1.thetachar
        ΔCp = db_filtered_1.DeltaCp
        φ₀ = db_filtered_1.phi0
        if "Cnumber" in names(db)
			Annotation = db_filtered_1.Annotation 
            #Dag = Array{Substance}(undef, length(Name))
            #for i=1:length(Name)
            #    Dag[i] = diffusivity(db_filtered_1.Molmass[i], 
            #                                    db_filtered_1.Cnumber[i], 
            #                                    db_filtered_1.Hnumber[i], 
            #                                    db_filtered_1.Onumber[i], 
            #                                    db_filtered_1.Nnumber[i], 
            #                                    db_filtered_1.Ringnumber[i], 
            #                                    gas)
            #end
        else # newer database version without information about the solute structure
            nCat = length(names(db))-8
            if nCat < 1
                Annotation = string.(db_filtered_1.Source)
            else
                Annotation = string.(db_filtered_1.Source)
                for i=1:nCat
                    for j=1:length(Annotation)
                        if typeof(db_filtered_1[j,8+i]) != Missing
                            Annotation[j] = string(Annotation[j], ", ", db_filtered_1[j,8+i])
                        end
                    end
                end
            end
            #Dag = Array{Substance}(undef, length(Name))
            #for i=1:length(Name)
            #    Dag[i] = diffusivity(CAS[i], gas)
            #end
        end
        # correct assignment of the t₀ and τ₀ values to the correct values from the input
        indices = Array{Int}(undef, length(Name))
        for i=1:length(Name)
            indices[i] = findfirst(db_filtered_1.CAS[i].==id.CAS)
        end
        t₀_ = t₀[indices]
        τ₀_ = τ₀[indices]  
	end
    sub = Array{Substance}(undef, length(Name))
    for i=1:length(Name)
        Dag = diffusivity(CAS[i], gas)
        sub[i] = Substance(Name[i],
                            CAS[i],
                            Tchar[i], 
                            θchar[i], 
                            ΔCp[i], 
                            φ₀[i],
                            Annotation[i],
                            Dag, 
                            t₀_[i],
                            τ₀_[i])
    end
	return sub
end

"""
    load_solute_database(db_path, db, sp, gas, solutes, t₀, τ₀) 

Load the data of `solutes` for the stationary phase `sp` and the mobile
phase `gas` from the database file `db` (located in `db_path`) into an array
of the structure `Substance`. 

# Arguments
* `db_path::String`: Path to the database file.
* `db::String`: Name of the database file. 
* `sp::String`: Name of the stationary phase.
* `gas::String`: Name of the mobile phase.
* `solutes::Array{<:AbstractString,1}`: Name of the solutes.
* `t₀::Array{Float64,1}`: Initial start times of the solutes.
* `τ₀::Array{Float64,1}`: Initial peak widths of the solutes. 

# Examples
```julia
julia> sub = load_solute_database("path/to/the/file", "db.csv", "DB5", "He", ["C10", "C11"], [0.0, 0.0], [0.5, 0.5])
```
"""
function load_solute_database(db_path::String, db::String, sp::String, gas::String, solutes::Array{<:AbstractString,1}, t₀::Array{Float64,1}, τ₀::Array{Float64,1})
	# load the information about a solute from a data base and collecting these informations in the 
    # structure Substance
    # 
	# db_path ... Path to the database file
	# db ... Name of the database
	# sp ... Name of the stationary phase
	# gas ... Name of the used mobile phase gas
	# solutes ... Names of the solutes for which the informations should be used
    # τ₀ ... initial values of the peak width for the solutes
    # t₀ ... initial values of the time for the solutes
	db_dataframe = DataFrame(CSV.File(string(db_path,"/",db), header=1, silencewarnings=true))
    sub = load_solute_database(db_dataframe, sp, gas, solutes, t₀, τ₀)
	return sub
end

"""
    all_solutes(sp, db) 

Extract the name of all solutes for which data in a database `db` and the
stationay phase `sp` is available. 

# Arguments
* `sp`: Name of the stationary phase.
* `db`: DataFrame of the database.

# Examples
```julia
julia> all = all_solutes("DB5", db)
```
"""
function all_solutes(sp, db)
	db_filter = filter([:Phase] => x -> x==sp, db)
	solutes = string.(db_filter.Name)
	return solutes
end

"""
    conventional_program(CP; time_unit="min")

Translate the conventional program notation into a vector of time steps and value steps (temperature, pressure, flow) used in GasChromatographySimulator.Program

The conventional temperature program is defined as an array of the following form (for a temperature program):
`CP = [T₁, t₁, r₁, T₂, t₂, r₂, T₃, t₃, ...]` corresponding to the notation:
`T₁(t₁) - r₁ - T₂(t₂) - r₂ - T₃(t₃) - ...` which can be read as:
Starting temperature `T₁` is holded for time `t₁`. After the holding time the temperature increases linearly with the heating rate `r₁`, until temperature `T₂` is reached. This temperature is held for the time `t₂` after which the temperature increases linearly by the heating rate `r₂` until temperature `T₃` is reached, which is hold for the time `t₃`, and so on. 

The option `time_unit` determines the unit of time in the program `CP`. For `time_unit = "min"` (default) the times are measured in minutes and the heating rates in °C/min. For `time_unit = "s"` the times are measured in seconds and the heating rate in °C/s. 
"""
function conventional_program(CP; time_unit="min")
    if time_unit == "min"
        c = 60.0
    elseif time_unit == "s"
        c = 1.0
    end
    value_steps = Float64[]
    for i=1:3:length(CP) # every third CP-entry starting from 1 is a value step
        push!(value_steps, CP[i])
        push!(value_steps, CP[i])
    end
    hold_times = Float64[]
    for i=2:3:length(CP)
        push!(hold_times, CP[i])
    end
    heating_rates = Float64[]
    for i=3:3:length(CP)
        push!(heating_rates, CP[i])
    end
    time_steps = Array{Float64}(undef, length(value_steps))
    time_steps[1] = 0.0
    for i=2:2:length(value_steps)
        time_steps[i] = hold_times[Int(i/2)] * c
        if i<length(value_steps)
            time_steps[i+1] = (value_steps[i+1] - value_steps[i])/heating_rates[Int(i/2)] * c
        end
    end
    return time_steps, value_steps
end

"""
    temperature_program(time_steps, value_steps; time_unit="min")

Translate the vector of time steps and value steps (temperature, pressure, flow) into a conventional program notation.

The conventional temperature program is defined as an array of the following form (for a temperature program):
`CP = [T₁, t₁, r₁, T₂, t₂, r₂, T₃, t₃, ...]` corresponding to the notation:
`T₁(t₁) - r₁ - T₂(t₂) - r₂ - T₃(t₃) - ...` which can be read as:
Starting temperature `T₁` is holded for time `t₁`. After the holding time the temperature increases linearly with the heating rate `r₁`, until temperature `T₂` is reached. This temperature is held for the time `t₂` after which the temperature increases linearly by the heating rate `r₂` until temperature `T₃` is reached, which is hold for the time `t₃`, and so on. 

The option `time_unit` determines the unit of time in the program `CP`. For `time_unit = "min"` (default) the times are measured in minutes and the heating rates in °C/min. For `time_unit = "s"` the times are measured in seconds and the heating rate in °C/s. 
"""
function temperature_program(time_steps, value_steps; time_unit="min")
    if time_unit == "min"
        c = 60.0
    elseif time_unit == "s"
        c = 1.0
    end
    # identify temperature pairs (the same value of temperature at neigboring elements)
    index_pair = Int[]
    for i=1:(length(value_steps)-1)
        if value_steps[i] == value_steps[i+1]
            push!(index_pair, i)
        end
    end
    # identify single temperatures
    index_single = Int[]
    for i=1:length(value_steps)
        if (i in index_pair) == false && (i in index_pair.+1) == false
            push!(index_single, i)
        end
    end
    values = value_steps[sort([index_pair; index_single])]
    # every (1+(i-1)*3)th element of VP is a value element of `values`
    
    thold = Array{Float64}(undef, length(values))
    a = sort([index_pair.+1; index_single])
    for i=1:length(values)
        if a[i] in index_single # holding time is zero for single entrys
            thold[i] = 0.0
        else # holding times are for paired temperatures the time_steps[index_pair.+1]
            thold[i] = time_steps[a[i]] / c
        end
    end 
    # every (2+(i-1)*3)th element of VP is a holding time

    rate = Array{Float64}(undef, length(values)-1)
    theat = time_steps[sort([index_pair; index_single])[2:end]]
    for i=1:(length(values)-1)
        rate[i] = (values[i+1] - values[i])/theat[i] * c
    end

    VP = Array{Float64}(undef, 2 + (length(values)-1)*3)
    for i=1:length(values)
        VP[1+(i-1)*3] = values[i]
    end
    for i=1:length(thold)
        VP[2+(i-1)*3] = thold[i]
    end
    for i=1:length(rate)
        VP[3+(i-1)*3] = rate[i]
    end
    return VP
end

"""
    common_time_steps(time_steps_1, time_steps_2)

Estimate a new set of time steps, which represents the combination of `time_steps_1` and `time_steps_2`.
"""
function common_time_steps(time_steps_1, time_steps_2)
	# constructs the new timesteps common for all moduls
	ctselements = Float64[]
    for j=1:length(time_steps_1)
        push!(ctselements, cumsum(time_steps_1)[j])
    end
    for j=1:length(time_steps_2)
        push!(ctselements, cumsum(time_steps_2)[j])
    end
	new_time_steps = [0.0; diff(sort(unique(ctselements)))]
	return new_time_steps
end

"""
    new_value_steps(value_steps, time_steps, new_time_steps)

Estimate the new value steps at the `new_time_steps` from the original set of `value_steps` over `time_steps`. The new values at new time steps are calculated from a linear change of the value between the original time steps.
""" 
function new_value_steps(value_steps, time_steps, new_time_steps)
    new_values = Array{Float64}(undef, length(new_time_steps))
    index_calc = Int[]
    index_old = Int[]
    for i=1:length(new_time_steps)
        if cumsum(new_time_steps)[i] in cumsum(time_steps)
            j = findfirst(cumsum(new_time_steps)[i].==cumsum(time_steps))
            new_values[i] = value_steps[j]
            push!(index_old, i)
        else
            push!(index_calc, i)
        end
    end
    for i=1:length(index_calc)
        i1 = findlast(index_old.<index_calc[i])
        i2 = findfirst(index_old.>index_calc[i])
        if isnothing(i1) || isnothing(i2)
            new_values[index_calc[i]] = value_steps[end]
        else
            if time_steps[i2] == 0.0
                i2 = i2 + 1
            end
            v1 = value_steps[i1]
            v2 = value_steps[i2]
            t1 = cumsum(time_steps)[i1]
            t2 = cumsum(time_steps)[i2]
            rate = (v2 - v1)/(t2 - t1)
            if t1 == t2
                new_values[index_calc[i]] = v2
            else
                new_values[index_calc[i]] = v1 + rate * (cumsum(new_time_steps)[index_calc[i]] - t1)
            end
        end
    end
    return new_values
end

#---End-Functions-used-for-Parameter-construction--- 

#---Begin-Solving-Functions---
"""
    simulate(par::Parameters)

Simulate the GC system defined by the structure `par`.

Note: Based on the option for `odesys` the result is different. For `odesys =
true` the result is a dataframe (the peaklist) and the solution of the ODEs
as a system (solution structure from DifferentialEquations.jl). If `odesys =
false` the result is a dataframe (the peaklist) and the two solutions of the
ODEs for time ``t(z)`` and peak variance ``τ²(z)``.
"""
function simulate(par)
    if par.opt.odesys==true
        sol = solve_system_multithreads(par)
    	pl = GasChromatographySimulator.peaklist(sol, par)
        return pl, sol
	else
		sol, peak = solve_multithreads(par)
    	pl = GasChromatographySimulator.peaklist(sol, peak, par)
        return pl, sol, peak
	end
end

"""
    solve_system_multithreads(par::Parameters)

Simulate the GC system defined by the structure `par` by solving the
ODEs for ``t(z)`` and ``τ²(z)`` together as a system of ODEs using multiple
threads (parallel computing) for the simulation of different solutes. 

Note: The result is an array of the solution structure from DifferentialEquations.jl.

# Examples

```julia
julia> sol = solve_system_multithreads(par)
```
"""
function solve_system_multithreads(par)
	n = length(par.sub)
	sol = Array{Any}(undef, n)
	Threads.@threads for i=1:n
		sol[i] = solving_odesystem_r(par.col, par.prog, par.sub[i], par.opt)
	end
	return sol
end

"""
    solve_multithreads(par::Parameters)

Simulate the GC system defined by the structure `par` by solving the
ODEs for ``t(z)`` and ``τ²(z)`` separatly (solving ``t(z)`` and using this result
to solve for ``τ²(z)``) using multiple threads (parallel computing) for the
simulation of different solutes.

Note: The result are two arrays of the solution structure from
DifferentialEquations.jl.

# Examples

```julia
julia> sol, peak = solve_multithreads(par)
```
"""
function solve_multithreads(par)
    n = length(par.sub)
    sol = Array{Any}(undef, n)
    peak = Array{Any}(undef, n)
    Threads.@threads for i=1:n
        sol[i] = solving_migration(par.col, par.prog, par.sub[i], par.opt)
        peak[i] = solving_peakvariance(sol[i], par.col, par.prog, par.sub[i], par.opt)
    end
    return sol, peak
end

"""
    solving_migration(col::Column, prog::Program, sub::Substance, opt::Options)

Solve for the migration ``t(z)`` of solute `sub` in the GC Column `col` with
the program `prog` and the options `opt`.

Note: The result is the solution structure from
DifferentialEquations.jl.
"""
function solving_migration(col::Column, prog::Program, sub::Substance, opt::Options)
	f_tz(t,p,z) = residency(z, t, prog.T_itp, prog.Fpin_itp, prog.pout_itp, col.L, col.d, col.df, col.gas, sub.Tchar, sub.θchar, sub.ΔCp, sub.φ₀; ng=opt.ng, vis=opt.vis, control=opt.control)
    t₀ = sub.t₀
    zspan = (0.0,col.L)
    prob_tz = ODEProblem(f_tz, t₀, zspan)
    solution_tz = solve(prob_tz, alg=opt.alg, abstol=opt.abstol,reltol=opt.reltol)
    return solution_tz
end

"""
    solving_peakvariance(solution_tz, col::Column, prog::Program, sub::Substance, opt::Options)

Solve for the development of the peak variance ``τ²(z)`` of solute `sub` in the GC Column `col` with
the program `prog` and the options `opt` during its migration defined by `solution_tz`.

Note: The result is the solution structure from
DifferentialEquations.jl.
"""
function solving_peakvariance(solution_tz, col, prog, sub, opt)
    t(z) = solution_tz(z)
    p = [col, prog, sub, opt]
    f_τ²z(τ²,p,z) = peakode(z, t(z), τ², col, prog, sub, opt)
    τ²₀ = sub.τ₀^2
    zspan = (0.0, col.L)
    prob_τ²z = ODEProblem(f_τ²z, τ²₀, zspan, p)
    solution_τ²z = solve(prob_τ²z, alg=opt.alg, abstol=opt.abstol,reltol=opt.reltol)
    return solution_τ²z
end

"""
    solving_odesystem_r(col::Column, prog::Program, sub::Substance, opt::Options)

Solve the migration ``t(z)`` and peak variance development ``τ²(z)`` of solute `sub` in the GC Column `col` with
the program `prog` and the options `opt` as a system of ODEs.

Note: The result is the solution structure from
DifferentialEquations.jl.

See also: [`odesystem_r!`](@ref)
"""
function solving_odesystem_r(col::Column, prog::Program, sub::Substance, opt::Options)
    t₀ = [sub.t₀; sub.τ₀^2]
    zspan = (0.0,col.L)
	p = [col, prog, sub, opt]
    prob = ODEProblem(odesystem_r!, t₀, zspan, p)

    solution = solve(prob, alg=opt.alg, abstol=opt.abstol,reltol=opt.reltol)

    if solution.t[end]<col.L
        solution = solve(prob, alg=opt.alg, abstol=opt.abstol,reltol=opt.reltol, dt=col.L/1000000)
    end
    return solution
end

"""
    odesystem_r!(dt, t, p, z)

The ODE system for migration ``t(z)`` and peak variance development
``τ²(z)``.

``\\frac{dt}{dz} = r(z, t(z))``

``\\frac{dτ^2}{dz} = H(z, t(z)) r(z, t(z)) + 2 τ^2(z, t(z))
\\frac{∂r}{∂t}(z,t(z))``

See also: [`solving_odesystem_r`](@ref), [`peakode`](@ref)
"""
function odesystem_r!(dt, t, p, z)
    # t[1] ... t time
    # t[2] ... τ² band variance
	col = p[1]
	prog = p[2]
	sub = p[3]
	opt = p[4]
    dt[1] = residency(z, t[1], prog.T_itp, prog.Fpin_itp, prog.pout_itp, col.L, col.d, col.df, col.gas, sub.Tchar, sub.θchar, sub.ΔCp, sub.φ₀; ng=opt.ng, vis=opt.vis, control=opt.control)
    dt[2] = peakode(z, t[1], t[2], col, prog, sub, opt)
end

"""
    peakode(z, t, τ², col, prog, sub, opt)

The second ODE function for the ODE system describing the peak variance
development ``τ²(z)``, using (in parts) automatic differentiation.

``\\frac{dτ^2}{dz} = H(z, t(z)) r(z, t(z)) + 2 τ^2(z, t(z))
\\frac{∂r}{∂t}(z,t(z))``

**TODO**: alternative to QuadGK.jl for integration which is available for
ForwardDiff.jl 

See also: [`solving_odesystem_r`](@ref), [`odesystem_r!`](@ref)
"""
function peakode(z, t, τ², col, prog, sub, opt)
	# alternative function
    if opt.ng==true
        r_ng(zt) = residency(zt[1], zt[2], prog.T_itp, prog.Fpin_itp, prog.pout_itp, col.L, col.d, col.df, col.gas, sub.Tchar, sub.θchar, sub.ΔCp, sub.φ₀; ng=true, vis=opt.vis, control=opt.control)
        H_ng(z,t) = plate_height(z, t, prog.T_itp, prog.Fpin_itp, prog.pout_itp, col.L, col.d, col.df, col.gas, sub.Tchar, sub.θchar, sub.ΔCp, sub.φ₀, sub.Dag; ng=true, vis=opt.vis, control=opt.control)
        ∂r∂t_ng(z,t) = ForwardDiff.gradient(r_ng, [z, t])[2]
        return H_ng(z,t)*r_ng([z,t])^2 + 2*τ²*∂r∂t_ng(z,t)
    else
		d(z) = col.d(z)
		T(z,t) = prog.T_itp(z, t)
		η(z,t) = GasChromatographySimulator.viscosity(z, t, prog.T_itp, col.gas; vis=opt.vis)
		c(zt) = η(zt[1], zt[2])*T(zt[1], zt[2])
		∂c∂t(z,t) = ForwardDiff.gradient(c, [z,t])[2]
        #if opt.control == "Pressure"
		#    pi2(t) = prog.Fpin_itp(t)^2
        #elseif opt.control == "Flow"
            pi2(t) = inlet_pressure(t, prog.T_itp, prog.Fpin_itp, prog.pout_itp, col.L, col.d, col.gas; vis=opt.vis, control=opt.control)^2
        #end
		po2(t) = prog.pout_itp(t)^2
		∂pi2∂t(t) = ForwardDiff.derivative(pi2, t)
		∂po2∂t(t) = ForwardDiff.derivative(po2, t)
        κ(z,t) = flow_restriction(z, t, prog.T_itp, col.d, col.gas; vis=opt.vis)
		κL(t) = κ(col.L,t)
		∂κ∂t(z,t) = quadgk(y -> d(y)^-4*∂c∂t(y,t), 0.0, z)[1]
        ∂κL∂t(t) = ∂κ∂t(col.L,t)
		e(z,t) = κ(z,t)/κL(t)
		∂e∂t(z,t) = (∂κ∂t(z,t)*κL(t)-κ(z,t)*∂κL∂t(t))/κL(t)^2
		p(z,t) = pressure(z, t, prog.T_itp, prog.Fpin_itp, prog.pout_itp, col.L, col.d, col.gas; vis=opt.vis, control=opt.control)
		∂p∂t(z,t) = 1/(2*p(z,t))*(∂pi2∂t(t)-(∂e∂t(z,t)*(pi2(t)-po2(t))+e(z,t)*(∂pi2∂t(t)-∂po2∂t(t))))
		rM(z,t) = mobile_phase_residency(z, t, prog.T_itp, prog.Fpin_itp, prog.pout_itp, col.L, col.d, col.gas; vis=opt.vis, control=opt.control)
		a(z,t) = κL(t)*p(z,t)
		∂a∂t(z,t) = ∂κL∂t(t)*p(z,t)+κL(t)*∂p∂t(z,t)
		b(zt) = T(zt[1],zt[2])*(pi2(zt[2])-po2(zt[2]))
		∂b∂t(z,t) = ForwardDiff.gradient(b, [z,t])[2]
		∂rM∂t(z,t) = 64*d(z)^2*((∂a∂t(z,t)*b([z,t])-a(z,t)*∂b∂t(z,t))/b([z,t])^2)
		
		k(zt) = retention_factor(zt[1], zt[2], prog.T_itp, col.d, col.df, sub.Tchar, sub.θchar, sub.ΔCp, sub.φ₀)
        ∂k∂t(z,t) = ForwardDiff.gradient(k, [z, t])[2]
		
		r(z,t) = residency(z, t, prog.T_itp, prog.Fpin_itp, prog.pout_itp, col.L, col.d, col.df, col.gas, sub.Tchar, sub.θchar, sub.ΔCp, sub.φ₀; vis=opt.vis, control=opt.control)
        H(z,t) = plate_height(z, t, prog.T_itp, prog.Fpin_itp, prog.pout_itp, col.L, col.d, col.df, col.gas, sub.Tchar, sub.θchar, sub.ΔCp, sub.φ₀, sub.Dag; vis=opt.vis, control=opt.control)
        ∂r∂t(z,t) = ∂rM∂t(z,t)*(1+k([z,t])) + rM(z,t)*∂k∂t(z,t)
        # AutoDiff problems with the intgral in flow_restriction κ and κL makes 
        # it nessecary to partly calculate the derivative manually
        # package Quadrature.jl makes it possible to differentiate for parameters
        # inside the integral, but not for the differentiation for upper or lower bound 
        # of the integral (issue #56)
        # if this can be solved, it should be straigthforward with AutoDiff
        return H(z,t)*r(z,t)^2 + 2*τ²*∂r∂t(z,t)
    end
end
#---End-Solving-Functions---

#---Begin-Result-Functions---
"""
    peaklist(sol, par)

Construct a DataFrame with the peak list of the solution `sol` of the
simulation of the GC system defined by `par`. 

# Output

The peaklist DataFrame consists of the entrys: 
* `Name`: Name of the solute.
* `tR`: Retention time of the solute (in s).
* `τR`: Peak width of the solute (in s). 
* `TR`: Temperature of the end of the column at the retention time (in °C).
* `σR`: Band width of the solute at retention time (in m).
* `uR`: Solute velocity at retention time (in m/s).
* `kR`: Retention factor of the solute at retention time.

# Examples

```julia
julia> pl = peaklist(sol, par)
...
```    
"""
function peaklist(sol, par)
	n = length(par.sub)
    # sol is solution from ODE system
    Name = Array{String}(undef, n)
    tR = Array{Float64}(undef, n)
    TR = Array{Float64}(undef, n)
    σR = Array{Float64}(undef, n)
    uR = Array{Float64}(undef, n)
    τR = Array{Float64}(undef, n)
    kR = Array{Float64}(undef, n)
    Res = fill(NaN, n)
    Threads.@threads for i=1:n
        Name[i] = par.sub[i].name
        if sol[i].t[end]==par.col.L
            tR[i] = sol[i].u[end][1]
            TR[i] = par.prog.T_itp(par.col.L, tR[i]) - 273.15 
            uR[i] = 1/residency(par.col.L, tR[i], par.prog.T_itp, par.prog.Fpin_itp, par.prog.pout_itp, par.col.L, par.col.d, par.col.df, par.col.gas, par.sub[i].Tchar, par.sub[i].θchar, par.sub[i].ΔCp, par.sub[i].φ₀; vis=par.opt.vis, control=par.opt.control)
            τR[i] = sqrt(sol[i].u[end][2])
            σR[i] = τR[i]*uR[i]
            kR[i] = retention_factor(par.col.L, tR[i], par.prog.T_itp, par.col.d, par.col.df, par.sub[i].Tchar, par.sub[i].θchar, par.sub[i].ΔCp, par.sub[i].φ₀)
        else
            tR[i] = NaN
            TR[i] = NaN
            uR[i] = NaN
            τR[i] = NaN
            σR[i] = NaN
            kR[i] = NaN
        end
    end  
    df = sort!(DataFrame(Name = Name, tR = tR, τR = τR, TR=TR, σR = σR, uR = uR, kR = kR, ), [:tR])
    Threads.@threads for i=1:n-1
        Res[i] = (df.tR[i+1] - df.tR[i])/(2*(df.τR[i+1] + df.τR[i]))
    end
    df[!, :Res] = Res  
    return df
end

"""
    peaklist(sol, peak, par)

Construct a DataFrame with the peak list of the solution `sol` and `peak` of
the simulation of the GC system defined by `par`. 

# Output

The peaklist DataFrame consists of the entrys: 
* `Name`: Name of the solute.
* `tR`: Retention time of the solute (in s).
* `τR`: Peak width of the solute (in s). 
* `TR`: Temperature of the end of the column at the retention time (in °C).
* `σR`: Band width of the solute at retention time (in m).
* `uR`: Solute velocity at retention time (in m/s).
* `kR`: Retention factor of the solute at retention time.

# Examples

```julia
julia> pl = peaklist(sol, peak, par)
...
```    
"""
function peaklist(sol, peak, par)
	n = length(par.sub)
    Name = Array{String}(undef, n)
    tR = Array{Float64}(undef, n)
    TR = Array{Float64}(undef, n)
    σR = Array{Float64}(undef, n)
    uR = Array{Float64}(undef, n)
    τR = Array{Float64}(undef, n)
    kR = Array{Float64}(undef, n)
    # add resolution
    Res = fill(NaN, n)
    Threads.@threads for i=1:n
        Name[i] = par.sub[i].name
        if sol[i].t[end]==par.col.L
            tR[i] = sol[i].u[end]
            TR[i] = par.prog.T_itp(par.col.L, tR[i]) - 273.15 
            uR[i] = 1/residency(par.col.L, tR[i], par.prog.T_itp, par.prog.Fpin_itp, par.prog.pout_itp, par.col.L, par.col.d, par.col.df, par.col.gas, par.sub[i].Tchar, par.sub[i].θchar, par.sub[i].ΔCp, par.sub[i].φ₀; vis=par.opt.vis, control=par.opt.control)
            τR[i] = sqrt(peak[i].u[end])
            σR[i] = τR[i]*uR[i]
            kR[i] = retention_factor(par.col.L, tR[i], par.prog.T_itp, par.col.d, par.col.df, par.sub[i].Tchar, par.sub[i].θchar, par.sub[i].ΔCp, par.sub[i].φ₀)
        else
            tR[i] = NaN
            TR[i] = NaN
            uR[i] = NaN
            τR[i] = NaN
            σR[i] = NaN
            kR[i] = NaN
        end
    end  
    df = sort!(DataFrame(Name = Name, tR = tR, τR = τR, TR=TR, σR = σR, uR = uR, kR = kR, ), [:tR])
    Threads.@threads for i=1:n-1
        Res[i] = (df.tR[i+1] - df.tR[i])/(2*(df.τR[i+1] + df.τR[i]))
    end
    df[!, :Res] = Res  
    return df
end

"""
    sol_extraction(sol, par)

Extract the points z=t, t=u1, τ²=u2 from the solution `sol` of
the ODE system of the GC system defined by `par` and exports them in a DataFrame.

# Examples

```julia
df_sol = sol_extraction(sol, par)
...
```    
"""
function sol_extraction(sol, par)
    # extract the points z=t, t=u1, τ²=u2 from the solution of
    # the ODE system
	n = length(par.sub)
    sol_z = Array{Any}(undef, n)
    sol_t = Array{Any}(undef, n)
    sol_τ² = Array{Any}(undef, n)
    solutes = Array{String}(undef, n)
    for i=1:n
        sol_z[i] = sol[i].t
        temp_t = Array{Float64}(undef, length(sol[i].t))
        temp_τ² = Array{Float64}(undef, length(sol[i].t))
        for j=1:length(sol[i].t)
                temp_t[j] = sol[i].u[j][1]
                temp_τ²[j] = sol[i].u[j][2]
        end
        sol_t[i] = temp_t
        sol_τ²[i] = temp_τ²
        solutes[i] = par.sub[i].name
    end
    df_sol = DataFrame(name=solutes, z=sol_z, t=sol_t, τ²=sol_τ²)
    return df_sol
end

"""
    sol_extraction(sol, peak, par)

Extract the points z_t=sol.t, t=sol.u, z_τ²=peak.t and τ²=peak.u from the
solution `sol` and `peak` of the ODEs of the GC system defined by `par` and exports them in a DataFrame.

# Examples

```julia
df_sol = sol_extraction(sol, peak, par)
...
```    
"""
function sol_extraction(sol, peak, par)
    # extract the points z=t, t=u, from the solution of
    # the first ODE (sol_tz)
    # and the points z=t, τ²=u, from the solution of 
    # the second ODE (peak_τz)
	n = length(par.sub)
    sol_z = Array{Any}(undef, n)
    sol_t = Array{Any}(undef, n)
    peak_z = Array{Any}(undef, n)
    peak_τ² = Array{Any}(undef, n)
    solutes = Array{String}(undef, n)
    for i=1:n
        sol_z[i] = sol[i].t
        sol_t[i] = sol[i].u
        peak_z[i] = peak[i].t
        peak_τ²[i] = peak[i].u
		solutes[i] = par.sub[i].name
	end
    df_sol = DataFrame(name=solutes, z_t=sol_z, t=sol_t, z_τ²=peak_z, τ²=peak_τ²)
    return df_sol
end

#----begin-notebooks-functions----------------------------------------------------------------------
## Physical-model-functions
include("./Model.jl")
## UI-functions
include("./UI.jl")
## plot-functions
include("./Plot.jl")
## misc-functions
include("./Misc.jl")
##---end-notebooks-functions---------------------------------------------------------------------------

end # module
