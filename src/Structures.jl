# Structures and Constructor functions
# ---Begin-Structures---
"""
    Column(L, d, a_d, df, a_df, sp, gas)

Structure describing the GC Column. 

# Arguments

* `L`: Length of the capillary measured in m (meter)
* `d`: A function `d(x, a_d)` of `x`, the position along the capillary, describing the diameter in m (meter). Or a number for a constant value.
* `a_d`: Parameters of the diameter function. 
* `d_f`: A function `d_f(x, a_df)` of `x`, describing the film thickness in m (meter). Or a number for a constant value.
* `a_df`: Parameters of the film thickness function. 
* `sp`: The name of the stationary phase.
* `gas`: The name of the mobile phase. Allowed values: He, H2 or N2.
"""
struct Column
    L                       # column length in m
    d                       # column internal diameter in m as function of x or a number
    a_d::Array{<:Real,1}   # parameters of the diameters function d(x)
    df                      # column film thickness in m as function of x or a number
    a_df::Array{<:Real}    # parameters of the film thickness function df(x)
    sp::String              # stationary phase of the column
    gas::String             # gas of the mobile phase ["He", "H2", "N2"]
    Column(L, d, a_d, df, a_df, sp, gas) = (gas!="He" && gas!="H2" && gas!="N2") ? error("Wrong selection for 'gas'. Choose 'He', 'H2' or 'N2'.") : new(L, d, a_d, df, a_df, sp, gas)
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
  	T_itp                                   # interpolation function of T(x,t)
    Fpin_itp                                # interpolation function of Fpin(t)
    pout_itp                                # interpolation function of pout(t)
    # add inner constructor to check the lengths of the Arrays and of the result of gf
    Program(ts, Ts, Fpis, pos, gf, agf, T, Fpin, po) = (length(ts)!=length(Ts) || length(ts)!=length(gf(0.0)) || length(ts)!=length(Fpis) || length(ts)!=length(pos)) || length(ts)!=size(agf)[1] ? error("Mismatch between length(time_steps) = $(length(ts)), length(temp_steps) = $(length(Ts)), length(Fpin_steps) = $(length(Fpis)), length(pout_steps) = $(length(pos)), length(gf(0.0)) = $(length(gf(0.0))) and size(a_gf)[1] = $(size(agf)[1]).") : new{typeof(gf)}(ts, Ts, Fpis, pos, gf, agf, T, Fpin, po)
end

"""
    Substance(name, CAS, Tchar, θchar, ΔCp, φ₀, ann, Cag, t₀, τ₀)

Structure to describe the properties of a solute, which migrates through the GC Column. These datas are in most cases read from a database with the function `load_solute_database()`.

# Arguments
* `name`: Name of the solute. 
* `CAS`: CAS number of the solute.
* `Tchar`: Characterisic temperature (in K). One of the three distribution-centric thermodynamic parameters describing the retention of this solute on the given stationary phase.
* `θchar`: Characterisic parameters (in °C). One of the three distribution-centric thermodynamic parameters describing the retention of this solute on the given stationary phase.
* `ΔCp`: Change of the isobaric heat capacity moving from the mobile to the stationary phase (in J mol⁻¹ K⁻¹). One of the three distribution-centric thermodynamic parameters describing the retention of this solute on the given stationary phase.
* `φ₀`: Dimensionless film thickness (φ ≈ df/d) of the column for which the thermodynamic parameters (Tchar, θchar, ΔCp) were estimated.
* `ann`: Annotations. In most cases the source of the data is noted here.
* `Cag`: The diffusitivity constant of the solute `a` in the mobile phase `g` (in...). It is calculated by the function `diffusitivity()`.
* `t₀`: Initial time of the solute (in s) at the start of the simulation.
* `τ₀`: Initial peak width of the solute (in s) at the start of the simulation. 

See also: [`load_solute_database`](@ref)
"""
struct Substance
    name::String        # name of solute
    CAS::Union{Missing,String}         # CAS registry number
    Tchar               # characteristic temperature in K
    θchar               # characteristic thermal constant in °C
    ΔCp                 # 3rd parameter
    φ₀                  # dimless film thickness for which Tchar, θchar and ΔCp were estimated
    ann::String         # annotations, e.g. the source of the data from which Tchar, θchar and ΔCp were estimated
    Cag                 # diffusitivity constant of analyt in a gas, calculate from structure (or from measurements)
    t₀                  # initial time in s  	
    τ₀                  # initial peak width in s   
end

"""
    Options(alg, abstol, reltol, Tcontrol, odesys, ng, vis, control, k_th)

Structure describing some general options for the simulation. 

# Arguments
* `alg`: The algorithm used for the ODE solver. Recommended: `OwrenZen5()` (default), `Tsit5()`, `Vern9()`, `BS5()` and `DP5()`. Legacy/alternative: `OwrenZen3()`, `OwrenZen4()`.
* `abstol`: The absolute tolerance for the ODE solver. Recommended value 1e-6 to 1e-8 (use 1e-8 for high-accuracy solvers like `Vern9()`).
* `reltol`: The relative tolerance for the ODE solver. Recommended value 1e-3 to 1e-5 (use 1e-5 for high-accuracy solvers like `Vern9()`).
* `Tcontrol`: Option defining at which point of the column the temperature program is calculated. The options are `inlet` (x=0) and `outlet` (x=L).
* `odesys`: Combine the ODEs for migration and peak-width into a system of ODEs (`odesys = true`) or solve the two ODEs separately (`odesys = false`).
* `ng`: Option to calculate the simulation without a gradient (`ng = true`) or with a gradient (`ng = false`).
* `vis`: Used model of viscosity. `HP` is a second-order polynomial taken from the HP flow calculator. `Blumberg` is an emperical formula according to the book
    `Temperature-programmed Gas Chromatography` by Leonid M. Blumberg (2010, Wiley-VCH).
* `control`: Control of the "Flow" or of the "Pressure" (at column inlet) during the program.
* `k_th`: Threshold for the maximum of the retention factor. If the calculated retention factor is bigger than `k_th` than the retention factor is set to the value `k_th`.
    This is done to avoid to small step widths in the solver for highly retained soultes at the beginning of a GC program. 

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
    k_th                # threshold for the max. possible retention factor
    # TODO: add check for the correct values of the options: Tcontrol, vis, control
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
    #d_func(x) = gradient(x, [d])
    #df_func(x) = gradient(x, [df])
    col = Column(L, d, [d], df, [df], sp, gas)
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
function Program(TP::Array{<:Real, 1}, FpinP::Array{<:Real, 1}, poutP::Array{<:Real, 1}, ΔTP::Array{<:Real, 1}, x₀P::Array{<:Real, 1}, L₀P::Array{<:Real, 1}, αP::Array{<:Real, 1}, Tcontrol::String, L::Real; time_unit="min")
    # using as gradient function the exponential model 'gf_exp(x,a_gf,Tcontrol)
    ts = Array{Array{Real,1}}(undef, 7)
    ts[1], Ts = GasChromatographySimulator.conventional_program(TP; time_unit=time_unit)
    ts[2], Fps = GasChromatographySimulator.conventional_program(FpinP; time_unit=time_unit)
    ts[3], pouts = GasChromatographySimulator.conventional_program(poutP; time_unit=time_unit)
    ts[4], ΔTs = GasChromatographySimulator.conventional_program(ΔTP; time_unit=time_unit)
    ts[5], x₀s = GasChromatographySimulator.conventional_program(x₀P; time_unit=time_unit)
    ts[6], L₀s = GasChromatographySimulator.conventional_program(L₀P; time_unit=time_unit)
    ts[7], αs = GasChromatographySimulator.conventional_program(αP; time_unit=time_unit)
    time_steps = Real[]
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
    a_gf = [zeros(length(time_steps)) zeros(length(time_steps)) Measurements.value(L).*ones(length(time_steps)) zeros(length(time_steps))]
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
* `pout`: Outlet pressure, "vacuum" (default), "atmosphere" or the outlet pressure in Pa(a).
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
    # remove additional 0.0 which are not at the first position
    ts1_ = ts1[[1; findall(0.0.!=ts1)]]
	Ts_ = Ts[[1; findall(0.0.!=ts1)]]
	ts2_ = ts2[[1; findall(0.0.!=ts2)]]
	Fps_ = Fps[[1; findall(0.0.!=ts2)]]
    time_steps = GasChromatographySimulator.common_time_steps(ts1_, ts2_)
    temp_steps = GasChromatographySimulator.new_value_steps(Ts_, ts1_, time_steps)
    Fpin_steps = GasChromatographySimulator.new_value_steps(Fps_, ts2_, time_steps)
    if pout == "vacuum"
        pout_steps = zeros(length(time_steps))
    elseif isa(pout, Number)
        pout_steps = pout.*ones(length(time_steps)) 
    else
        pout_steps = pn.*ones(length(time_steps))
    end
    a_gf = [zeros(length(time_steps)) zeros(length(time_steps)) Measurements.value(L).*ones(length(time_steps)) zeros(length(time_steps))]
    gf(x) = GasChromatographySimulator.gradient(x, a_gf)
    T_itp = GasChromatographySimulator.temperature_interpolation(time_steps, temp_steps, gf, L)
    Fpin_itp = GasChromatographySimulator.steps_interpolation(time_steps, Fpin_steps)
    pout_itp = GasChromatographySimulator.steps_interpolation(time_steps, pout_steps)
    prog = GasChromatographySimulator.Program(time_steps, temp_steps, Fpin_steps, pout_steps, gf, a_gf, T_itp, Fpin_itp, pout_itp)
    return prog
end

"""
    Options(;alg=OwrenZen5(), abstol=1e-6, reltol=1e-3, Tcontrol="inlet", odesys=true, ng=false, vis="Blumberg", control="Pressure", k_th=1e12)

Construct the structure `Options` with default values. 

# Arguments
* `alg`: The algorithm used for the ODE solver. Recommended: `OwrenZen5()` (default), `Tsit5()`, `Vern9()`, `BS5()` and `DP5()`. Legacy/alternative: `OwrenZen3()`, `OwrenZen4()`.
* `abstol`: The absolute tolerance for the ODE solver. Recommended value 1e-6 to 1e-8 (use 1e-8 for high-accuracy solvers like `Vern9()`).
* `reltol`: The relative tolerance for the ODE solver. Recommended value 1e-3 to 1e-5 (use 1e-5 for high-accuracy solvers like `Vern9()`).
* `Tcontrol`: Option defining at which point of the column the temperature
    program is calculated. The options are `inlet` (x=0) and `outlet` (x=L).
* `odesys`: Combine the ODEs for migration and peak-width into a system of
    ODEs (`odesys = true`) or solve the two ODEs separately (`odesys = false`).
* `ng`: Option to calculate the simulation without a gradient (`ng = true`)
    or with a gradient (`ng = false`).
* `vis`: Used model of viscosity. `HP` is a second-order polynomial taken from the HP flow calculator. `Blumberg` is an emperical formula according to the book
    `Temperature-programmed Gas Chromatography` by Leonid M. Blumberg (2010, Wiley-VCH) 
* `control`: Control of the "Flow" or of the "Pressure" (at column inlet) during the program
* `k_th`: Threshold for the maxima allowed value of retention factor.

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
function Options(;alg=OwrenZen5(), abstol=1e-6, reltol=1e-3, Tcontrol="inlet", odesys=true, ng=false, vis="Blumberg", control="Pressure", k_th=1e12)
    opt = Options(alg, abstol, reltol, Tcontrol, odesys, ng, vis, control, k_th)
    return opt
end

"""
    Options(alg, abstol, reltol, Tcontrol, odesys; ng=false, vis="Blumberg", control="Pressure", k_th=1e12)

Construct the structure `Options` with given values. 

# Arguments
* `alg`: The algorithm used for the ODE solver. Recommended: `OwrenZen5()` (default), `Tsit5()`, `Vern9()`, `BS5()` and `DP5()`. Legacy/alternative: `OwrenZen3()`, `OwrenZen4()`.
* `abstol`: The absolute tolerance for the ODE solver. Recommended value 1e-6 to 1e-8 (use 1e-8 for high-accuracy solvers like `Vern9()`).
* `reltol`: The relative tolerance for the ODE solver. Recommended value 1e-3 to 1e-5 (use 1e-5 for high-accuracy solvers like `Vern9()`).
* `Tcontrol`: Option defining at which point of the column the temperature
    program is calculated. The options are `inlet` (x=0) and `outlet` (x=L).
* `odesys`: Combine the ODEs for migration and peak-width into a system of
    ODEs (`odesys = true`) or solve the two ODEs separately (`odesys = false`).
* `ng`: Option to calculate the simulation without a gradient (`ng = true`)
    or with a gradient (`ng = false`).
* `vis`: Used model of viscosity. `HP` is a second-order polynomial taken from the HP flow calculator. `Blumberg` is an emperical formula according to the book
    `Temperature-programmed Gas Chromatography` by Leonid M. Blumberg (2010, Wiley-VCH) 
* `control`: Control of the "Flow" or of the "Pressure" (at column inlet) during the program
* `k_th`: Threshold for the maxima allowed value of retention factor.

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
function Options(alg, abstol, reltol, Tcontrol, odesys; ng=false, vis="Blumberg", control="Pressure", k_th=1e12)
    opt = Options(alg, abstol, reltol, Tcontrol, odesys, ng, vis, control, k_th)
    return opt
end

# Aliases for compatibility
constructor_System = Column
constructor_Program = Program
constructor_Options = Options
#---End-Constructor-functions-for-Structures---