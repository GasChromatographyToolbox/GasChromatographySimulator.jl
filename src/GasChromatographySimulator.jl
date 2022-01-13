module GasChromatographySimulator

using Reexport
@reexport using Interpolations
using QuadGK
@reexport using DifferentialEquations
using ForwardDiff
@reexport using DataFrames
@reexport using CSV

# some constants
Tst = 273.15            # K
R = 8.31446261815324    # J mol⁻¹ K⁻¹
Tn = 25.0 + Tst         # K
pn = 101300             # Pa

# ---Begin-Structures---
"""
    System(L, d, a_d, df, a_df, sp, gas)

Structure describing the GC system. 

# Arguments

* `L`: Length of the capillary measured in m (meter)
* `d`: A function `d(x, a_d)` of `x`, the position along the capillary, describing the diameter in m (meter).
* `a_d`: Parameters of the diameter function. 
* `d_f`: A function `d_f(x, a_df)` of `x`, describing the film thickness in m (meter).
* `sp`: The name of the stationary phase.
* `gas`: The name of the mobile phase. Allowed values: He, H2 or N2.
"""
struct System{Fd<:Function, Fdf<:Function}
    L                       # column length in m
    d::Fd                   # column internal diameter in m as function of x
    a_d::Array{Float64,1}   # parameters of the diameters function d(x)
    df::Fdf                 # column film thickness in m as function of x
    a_df::Array{Float64}    # parameters of the film thickness function df(x)
    sp::String              # stationary phase of the column
    gas::String             # gas of the mobile phase ["He", "H2", "N2"]
    System(L, d, a_d, df, a_df, sp, gas) = (gas!="He" && gas!="H2" && gas!="N2") ? error("Wrong selection for 'gas'. Choose 'He', 'H2' or 'N2'.") : new{typeof(d), typeof(df)}(L, d, a_d, df, a_df, sp, gas)
end

"""
    Program(time_steps, temp_steps, pin_steps,
    pout_steps, gf, a_gf, T_itp, pin_itp, pout_itp)

Structure to describe the temperature and pressure program of a GC system. The
function `gf` describes an optional thermal gradient.

# Arguments
* `time_steps`: Time steps in s (seconds). 
* `temp_steps`: Temperature steps in °C (degree celsius).
* `pin_steps`: Inlet pressure steps in Pa(a) (pascal, absolute).
* `pout_steps`: Outlet pressure steps in Pa(a) (pascal, absolute).
* `gf`: Gradient function `gf(x, a_gf)`, describes the thermal gradient.
* `a_gf`: Parameters of the gradient function.
* `T_itp`: Interpolated (linear) temperature `T(x,t)`, constructed from `time_steps`, `temp_steps` and `gf`.
* `pin_itp`: Interpolated (linear) inlet pressure `pin(t)`, constructed from `time_steps` and `pin_steps`.
* `pout_itp`: Interpolated (linear) outlet pressure `pout(t)`, constructed from `time_steps` and `pout_steps`.  

Note: The length of the arrays `time_steps`, `temp_steps`, `pin_steps` and `a_gf`
have to be the same.
"""
struct Program{Fgf<:Function}
    time_steps::Array{<:Real, 1}            # vector time steps for the temperature program
    temp_steps::Array{<:Real, 1}            # vector temperature steps for the temperature program 
    pin_steps::Array{<:Real, 1}             # vector inlet pressure steps for the pressure program
    pout_steps::Array{<:Real, 1}            # vector outlet pressure steps for the pressure program
    gf::Fgf                                 # function of x of the gradient form
    a_gf::Array{<:Real}                     # parameters of the gradient function gf(x)
  	T_itp::Interpolations.Extrapolation     # interpolation function of T(x,t)
    pin_itp::Interpolations.Extrapolation   # interpolation function of pin(t)
    pout_itp::Interpolations.Extrapolation	# interpolation function of pout(t)
    # add inner constructor to check the lengths of the Arrays and of the result of gf
    Program(ts, Ts, pis, pos, gf, agf, T, pin, po) = (length(ts)!=length(Ts) || length(ts)!=length(gf(0.0)) || length(ts)!=length(pis) || length(ts)!=length(pos)) || length(ts)!=size(agf)[1] ? error("Mismatch between length(time_steps) = $(length(ts)), length(temp_steps) = $(length(Ts)), length(pin_steps) = $(length(pis)), length(pout_steps) = $(length(pos)), length(gf(0.0)) = $(length(gf(0.0))) and size(a_gf)[1] = $(size(agf)[1]).") : new{typeof(gf)}(ts, Ts, pis, pos, gf, agf, T, pin, po)
end

"""
    Substance(name, CAS, Tchar, θchar, ΔCp, φ₀, ann, Dag, t₀, τ₀)

Structure to describe the properties of a solute, which migrates through the
GC system. These datas are in most cases read from a database with the
function `load_solute_database()`.

# Arguments
* `name`: Name of the solute. 
* `CAS`: CAS number of the solute.
* `Tchar`: Characterisic temperature (in K). One of the three
  distribution-centric thermodynamic parameters describing the retention of this
  solute on the given stationary phase.
* `θchar`: Characterisic parameters (in °C). One of the three
distribution-centric thermodynamic parameters describing the retention of
this solute on the given stationary phase.
* `ΔCp`: Change of the isobaric heat capacity moving from the mobile to the
    stationary phase (in J mol⁻¹ K⁻¹). One of the three
    distribution-centric thermodynamic parameters describing the retention of
    this solute on the given stationary phase.
* `φ₀`: Dimensionless film thickness (φ ≈ df/d) of the column for which the
    thermodynamic parameters (Tchar, θchar, ΔCp) were estimated.
* `ann`: Annotations. In most cases the source of the data is noted here.
* `Dag`: The diffusitivity of the solute `a` in the mobile phase `g` (in
    ...). It is calculated by the function `diffusitivity()`.
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
    Options(alg, abstol, reltol, Tcontrol, odesys, ng)

Structure describing some general options for the simulation. 

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
    integrals, e.g. in the `flow_restriction()` function. -> **TODO**: test
    package Quadrature.jl as alternative to QuadGK.jl for integration)

**TODO**: add option for the retention model ('ABC', 'K-centric')

For more informations about the arguments `alg`, `abstol` and `reltol` see
the documentation of the DifferentialEquations.jl package.
"""
struct Options
    alg                 # algorithmen for the ODE solver
    abstol              # absolute tolerance for ODE solver
    reltol              # relative tolerance for ODE solver 
    Tcontrol::String    # temperature control at 'inlet' (top) or 'outlet' (bottom) of the column
	odesys::Bool  		# calculate the two ODEs (migration and peak-width) separately (false) or 
                        # combined as a system of ODEs (true)
                        # ??? add 'ng' as option ???
                        # for compatibility with previous code make constructor
                        # Options(alg, abstol, reltol, Tcontrol, odesys) =
                        # Options(alg, abstol, reltol, Tcontrol, odesys, false)
                        # (ng=false as default)
    ng::Bool
end

"""
    Parameters(sys, prog, sub, opt)

Structure describing all parameters for the simulation of a GC system. 

# Arguments
* `sys`: Structure `Systems` describing the parameters of the GC column and
    mobile phase gas.
* `prog`: Structure `Program` describing the temperature and pressure
    program of a GC system.
* `sub`: An array of the structure `Substance` describing the parameters of
    the solutes which are separated in the GC simulation. 
* `opt`: Structure `Options` describing additional option parameters.
"""
struct Parameters
    sys::System
    prog::Program
    sub::Array{Substance,1}
    opt::Options
end
# ---End-Structures---

# ---Begin-Constructor-functions-for-Structures---

"""
    System(L, d, df, sp, gas)

Construct the structure `Systems` with given values for the case
of constant diameter `d` and film thickness `df`. 

# Arguments
* `L`: Length of the capillary measured in m (meter).
* `d`: Diameter of the capillary measured in m (meter). 
* `d_f`: Film thickness of the capillary measured in m (meter).
* `sp`: The name of the stationary phase.
* `gas`: The name of the mobile phase. Allowed values: He, H2 or N2.

# Examples
```julia
julia> System(10.0, 0.1e-3, 0.1e-6, "DB5", "He")
	```
"""
function System(L, d, df, sp, gas)
    # function to construct the System structure
    # for the case of constant diameter and constant film thickness
    d_func(x) = gradient(x, [d])
    df_func(x) = gradient(x, [df])
    sys = System(L, d_func, [d], df_func, [df], sp, gas)
    return sys
end

"""
    Program(time_steps, temp_steps, pin_steps, pout_steps, a_gf, Tcontrol, L)

Construct the structure `Program` with given values. 

# Arguments
* `time_steps`: Time steps in s (seconds). 
* `temp_steps`: Temperature steps in °C (degree celsius).
* `pin_steps`: Inlet pressure steps in Pa(a) (pascal, absolute).
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
function Program(time_steps::Array{<:Real, 1}, temp_steps::Array{<:Real, 1}, pin_steps::Array{<:Real, 1}, pout_steps::Array{<:Real, 1}, a_gf::Array{<:Real, 2}, Tcontrol, L)
    # function to construct the Program structure
    # using as gradient function the exponential model 'gf_exp(x,a_gf,Tcontrol)'
    gf(x) = gradient(x, a_gf; Tcontrol=Tcontrol)
    T_itp = temperature_interpolation(time_steps, temp_steps, gf, L)
    pin_itp = pressure_interpolation(time_steps, pin_steps)
    pout_itp = pressure_interpolation(time_steps, pout_steps)
    prog = Program(time_steps, temp_steps, pin_steps, pout_steps, gf, a_gf, T_itp, pin_itp, pout_itp)
    return prog
end

"""
    Program(time_steps, temp_steps, pin_steps, pout_steps, ΔT_steps, x₀_steps, L₀_steps, α_steps, Tcontrol, L) 

Construct the structure `Program` with given values. 

# Arguments
* `time_steps`: Time steps in s (seconds). 
* `temp_steps`: Temperature steps in °C (degree celsius).
* `pin_steps`: Inlet pressure steps in Pa(a) (pascal, absolute).
* `pout_steps`: Outlet pressure steps in Pa(a) (pascal, absolute).
* `ΔT_steps`: Parameters of the gradient function. Temperature difference in °C.
* `x₀_steps`: Parameters of the gradient function. Spatial offset of the gradient in m.
* `L₀_steps`: Parameters of the gradient function. Distance over which the temperature difference in ΔT_steps is measured, in m.
* `α_steps`: Parameters of the gradient function. Factor describing the gradient profile. 
* `Tcontrol`: Option defining at which point of the column the temperature
program is calculated. The options are `inlet` (x=0) and `outlet` (x=L).
* `L`: Length of the capillary measured in m (meter).

The length of the arrays `time_steps`, `temp_steps`, `pin_steps`, `pout_steps`, `ΔT_steps`, `x₀_steps`, `L₀_steps` and `α_steps`
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
function Program(time_steps::Array{<:Real, 1}, temp_steps::Array{<:Real, 1}, pin_steps::Array{<:Real, 1}, pout_steps::Array{<:Real, 1}, ΔT_steps::Array{<:Real, 1}, x₀_steps::Array{<:Real, 1}, L₀_steps::Array{<:Real, 1}, α_steps::Array{<:Real, 1}, Tcontrol, L)
    # function to construct the Program structure
    # using as gradient function the exponential model 'gf_exp(x,a_gf,Tcontrol)'
    a_gf = [ΔT_steps x₀_steps L₀_steps α_steps]
    gf(x) = gradient(x, a_gf; Tcontrol=Tcontrol)
    T_itp = temperature_interpolation(time_steps, temp_steps, gf, L)
    pin_itp = pressure_interpolation(time_steps, pin_steps)
    pout_itp = pressure_interpolation(time_steps, pout_steps)
    prog = Program(time_steps, temp_steps, pin_steps, pout_steps, gf, a_gf, T_itp, pin_itp, pout_itp)
    return prog
end

"""
    Program(time_steps, temp_steps, pin_steps, pout_steps, L)

Construct the structure `Program` with given values for the case
without a thermal gradient. 

# Arguments
* `time_steps::Array{<:Real, 1}`: Time steps in s (seconds). 
* `temp_steps::Array{<:Real, 1}`: Temperature steps in °C (degree celsius).
* `pin_steps::Array{<:Real, 1}`: Inlet pressure steps in Pa(a) (pascal, absolute).
* `pout_steps::Array{<:Real, 1}`: Outlet pressure steps in Pa(a) (pascal, absolute).
* `L`: Length of the capillary measured in m (meter).

The length of the arrays `time_steps`, `temp_steps`, `pin_steps` and `pout_steps`
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
function Program(time_steps::Array{<:Real, 1}, temp_steps::Array{<:Real, 1}, pin_steps::Array{<:Real, 1}, pout_steps::Array{<:Real, 1}, L)
    # function to construct the Program structure
    # without a thermal gradient
    # using as gradient function the exponential model 'gf_exp(x,a_gf,Tcontrol)'
    a_gf = [zeros(length(time_steps)) zeros(length(time_steps)) L.*ones(length(time_steps)) zeros(length(time_steps))]
    gf(x) = gradient(x, a_gf)
    T_itp = temperature_interpolation(time_steps, temp_steps, gf, L)
    pin_itp = pressure_interpolation(time_steps, pin_steps)
    pout_itp = pressure_interpolation(time_steps, pout_steps)
    prog = Program(time_steps, temp_steps, pin_steps, pout_steps, gf, a_gf, T_itp, pin_itp, pout_itp)
    return prog
end

"""
    Options(;alg=OwrenZen5(), abstol=1e-6, reltol=1e-3, Tcontrol="inlet", odesys=true, ng=false)

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
function Options(;alg=OwrenZen5(), abstol=1e-6, reltol=1e-3, Tcontrol="inlet", odesys=true, ng=false)
    opt = Options(alg, abstol, reltol, Tcontrol, odesys, ng)
    return opt
end

"""
    Options(alg, abstol, reltol, Tcontrol, odesys; ng=false)

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

For more informations about the arguments `alg`, `abstol` and `reltol` see
the documentation of the DifferentialEquations.jl package.

# Examples
```julia
julia> Options(OwrenZen3(), 1e-7, 1e-4, "inlet", true)
```

```julia
julia> Options(OwrenZen3(), 1e-7, 1e-4, "inlet", true; ng=true)
```
"""
function Options(alg, abstol, reltol, Tcontrol, odesys; ng=false)
    opt = Options(alg, abstol, reltol, Tcontrol, odesys, ng)
    return opt
end

# Aliases for compatibility
constructor_System = System
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
    pressure_interpolation(time_steps, press_steps)

Construct the pressure function depending on time `t`.  

# Arguments
* `time_steps::Array{<:Real,1}`: Time steps in s (seconds). 
* `press_steps::Array{<:Real,1}`: Pressure steps in Pa (Pascal).

The pressure between the `time_steps` is linear interpolated between the
corresponding values of `press_steps`  

# Examples
```julia
julia> pin_itp = pressure_interpolation([0.0, 60.0, 300.0, 120.0], 
                                    [300000.0, 300000.0, 400000.0, 400000.0])
```
"""
function pressure_interpolation(time_steps::Array{<:Real,1}, press_steps::Array{<:Real,1})
    p_itp = LinearInterpolation((cumsum(time_steps), ), press_steps, extrapolation_bc=Flat())
    return p_itp
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
	# load the information about a solute from a data base and collecting these informations in the 
    # structure SubstanceGC2
    # 
	# db ... Dataframe of the database
	# sp ... Name of the stationary phase
	# gas ... Name of the used mobile phase gas
	# solutes ... Names of the solutes for which the informations should be used
    # τ₀ ... initial values of the peak width for the solutes
    # t₀ ... initial values of the time for the solutes
	#
	# at the moment only the new database config is supported (new data is appended in new rows)
	# also, if for a solute no data is available in the database no error is given!!!
	if size(db)[2]>14
		# old database config (additional stationary phases are appended in new columns, only one row for a solute)
        error("Data format not supported. Use the appended database structure.")
	else
		# new database config (simpler, additional stationary phases are appended in new rows, multiple entrys for a solute)
		# 1. Filter the stationary phase
		phases_db = unique(db.Phase)
		if isa(findfirst(phases_db.==sp), Nothing) && sp!=""
			error("Unknown selction of stationary phase. Choose one of these: $phases_db")
		elseif sp=="" # no stationary phase is selected, like for transferlines
			db_filtered = db[!,1:8] # use only the data of the first 8 columns
			# 2. Filter the solutes, not using multiple entrys
			db_filtered_1=unique(db_filtered[in(solutes).(db_filtered.Name),:])
			# use placeholder values
			Tchar = ones(length(solutes))
			θchar = ones(length(solutes))
			ΔCp = ones(length(solutes))
			φ₀ = ones(length(solutes))
			Annotation = fill("no sp", length(solutes))
		else
			db_filtered = filter([:Phase] => x -> x==sp, db)
			# 2. Filter the solutes
			db_filtered_1=db_filtered[in(solutes).(db_filtered.Name),:]
			# values
			Tchar = db_filtered_1.Tchar.+Tst
			θchar = db_filtered_1.thetachar
			ΔCp = db_filtered_1.DeltaCp
			φ₀ = db_filtered_1.phi0
			Annotation = db_filtered_1.Annotation
		end
        sub = Array{Substance}(undef, length(solutes))
        for i=1:length(solutes)
			Dag = diffusivity(db_filtered_1.Molmass[i], 
											db_filtered_1.Cnumber[i], 
											db_filtered_1.Hnumber[i], 
											db_filtered_1.Onumber[i], 
											db_filtered_1.Nnumber[i], 
											db_filtered_1.Ringnumber[i], 
											gas)
			sub[i] = Substance(db_filtered_1.Name[i],
										db_filtered_1.CAS[i],
										Tchar[i], 
										θchar[i], 
										ΔCp[i], 
										φ₀[i],
										Annotation[i],
										Dag, 
										t₀[i],
										τ₀[i])
		end
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
    diffusivity(M, Cn, Hn, On, Nn, Rn, gas)

Calculate the diffusivity `Dag` of solute `a` in gas `g` using the
emperical Fuller-Schettler-Giddings model [1].

[1] Fuller, Edward N.; Ensley, Keith; Giddings, J. Calvin, Diffusion of
Halogenated Hydrocarbons in Helium. The Effect of Structure on Collision
Cross Sections, The Journal of Physical Chemistry, Volume 73, Issue 11,
1969, 3679–3685

# Arguments
* `M`: Molar mass of the solute.
* `Cn`: Number of carbon atoms of the solute.
* `Hn`: Number of hydrogen atoms of the solute.
* `On`: Number of oxygen atoms of the solute.
* `Nn`: Number of nitrogen atoms of the solute.
* `Rn`: Number of closed rings of the structure of the solute.
* `gas`: The name of the mobile phase. Allowed values: He, H2 or N2.
"""
function diffusivity(M, Cn, Hn, On, Nn, Rn, gas)
    # calculates diffusitivity Dag of an analyte in a gas
    # from the (simplified) molecular formula of the solute
    # using the Fuller-Schettler-Giddings equation
    # 
    # M ... molar mass
    # Cn ... number of carbons
    # Hn ... number of hydrogens
    # On ... number of oxygens
    # Nn ... number of nitrogens
    # Rn ... number of ring structures
    # gas ... the used mobile phase gas
    # **TODO**: add other atoms (P, Cl, Si, ...)
    if gas=="H2"
        Vg = 6.12
        Mg = 2.02
    elseif gas=="He"
        Vg = 2.67
        Mg = 4
    elseif gas=="N2"
        Vg = 18.5
        Mg = 28.01
    elseif gas=="Ar"
        Vg = 16.2
        Mg = 39.95
    else
        error("Unknown selection of gas. Choose one of these: He, H2, N2 or Ar.")
    end
    Va = 15.9*Cn + 2.31*Hn + 6.11*On + 4.54*Nn - 18.3*Rn
    Dag = pn*sqrt(1/M+1/Mg)/(Vg^(1/3)+Va^(1/3))^2*1e-7 # m²/s
    return Dag
end
#---End-Functions-used-for-Parameter-construction--- 

#---Begin-Functions-of-the-physical-model---
"""
    pressure(x, t, T_itp, pin_itp, pout_itp, L, d, gas; ng=false)

Calculate the pressure at position `x` at time `t`.

# Arguments
* `x`: Position along the GC column, in m.
* `t`: Time in s.
* `T_itp`: Interpolated (linear) temperature `T(x,t)`.
* `pin_itp`: Interpolated (linear) inlet pressure `pin(t)`.
* `pout_itp`: Interpolated (linear) outlet pressure `pout(t)`.
* `L`: Length of the capillary measured in m (meter).
* `d`: Diameter of the GC column, in m. Can be a function of position `x`.
* `gas`: Name of the mobile phase gas.
* `ng`: Option to calculate the simulation without a gradient (`ng = true`,
    eq. 2)
or with a gradient (`ng = false`, eq. 1).

``p(x,t) =
\\sqrt(p_{in}(t)^2-\\frac{κ(x,t)}{κ_L(t)}\\left(p_{in}^2-p_{out}^2\\right))``
Eq. 1

``p(x,t) =
\\sqrt(p_{in}(t)^2-\\frac{x}{L}\\left(p_{in}^2-p_{out}^2\\right))`` Eq. 2

with ``κ(x,t)`` the flow restriction up to position `x` at time `t` and
``κ_L(t) = κ(x=L,t)`` the flow restriction of the whole column at
time `t`.

See also: [`flow_restriction`](@ref)
"""
function pressure(x, t, T_itp, pin_itp, pout_itp, L, d, gas; ng=false)
    if ng==true
        pp = sqrt(pin_itp(t)^2 - x/L*(pin_itp(t)^2-pout_itp(t)^2))
    else
        pp = sqrt(pin_itp(t)^2 - flow_restriction(x, t, T_itp, d, gas)/flow_restriction(L, t, T_itp, d, gas)*(pin_itp(t)^2-pout_itp(t)^2))
    end
    return pp
end

"""
    flow_restriction(x, t, T_itp, d, gas; ng=false)

Calculate the flow restriction ``κ`` up to position `x` at time `t`.

# Arguments
* `x`: Position along the GC column, in m.
* `t`: Time in s.
* `T_itp`: Interpolated (linear) temperature `T(x,t)`.
* `d`: Diameter of the GC column, in m. Can be a function of position `x`.
* `gas`: Name of the mobile phase gas.
* `ng`: Option to calculate the simulation without a gradient (`ng = true`,
    eq. 2)
or with a gradient (`ng = false`, eq. 1).

``κ(x,t) = \\int_0^x \\frac{η(y,t) T(y,t)}{d(y)^4}dy``
Eq. 1

``κ(x,t) = \\frac{η(t) T(t) x}{d^4}`` Eq. 2

with ``η(x,t)`` the viscosity of the mobile phase gas.

See also: [`viscosity`](@ref)
"""
function flow_restriction(x, t, T_itp, d, gas; ng=false)
    if ng==true
        κ = x*viscosity(x, t, T_itp, gas)*T_itp(x, t)*d(x)^-4
    else
        κ = quadgk(y -> viscosity(y, t, T_itp, gas)*T_itp(y, t)*d(y)^-4, 0, x)[1]
    end
    return κ
end

"""
    viscosity(x, t, T_itp, gas)

Calculate the (dynamic) viscosity of the mobile phase gas at position `x`
at time `t` in Pa s.

# Arguments
* `x`: Position along the GC column, in m.
* `t`: Time in s.
* `T_itp`: Interpolated (linear) temperature `T(x,t)`.
* `gas`: Name of the mobile phase gas.

``η(x,t) = η_{st}\\left(\\frac{T(x,t)}{T_{st}}\\right)^{(ξ_0 + ξ_1 \\frac{T(x,t)-T_{st}}{T_{st}})}`` 

with ``η_{st}``, ``ξ_0`` and ``ξ_1`` parameters dependent on the
mobile phase gas [1].

[1] Blumberg, Leonid M., Temperature-Programmed Gas Chromatography,
Wiley-VCH, 2010.
"""
function viscosity(x, t, T_itp, gas)
    # using empiric model from Blumberg.2010
    if gas=="He"
        ηst = 18.63e-6
        ξ₀ = 0.6958
        ξ₁ = -0.0071
    elseif gas=="H2"
        ηst = 8.382e-6
        ξ₀ = 0.6892
        ξ₁ = 0.005
    elseif gas=="N2"
        ηst = 16.62e-6
        ξ₀ = 0.7665
        ξ₁ = -0.0378
    elseif gas=="Ar"
        ηst = 21.04e-6
        ξ₀ = 0.8131
        ξ₁ = -0.0426
    else
        error("Unknown selection of gas. Choose one of these: He, H2, N2 or Ar.")
    end
    T = T_itp(x, t)
    η = ηst*(T/Tst)^(ξ₀ + ξ₁*(T-Tst)/Tst)
    return η
end

"""
    viscosity(T, gas)

Calculate the (dynamic) viscosity of the mobile phase gas at temperature `T` in Pa s.

# Arguments
* `T`: Temperature in K.
* `gas`: Name of the mobile phase gas.

``η(x,t) = η_{st}\\left(\\frac{T)}{T_{st}}\right)^{(ξ_0 + ξ_1 \\frac{T-T_{st}}{T_{st}})}`` 

with ``η_{st}``, ``ξ_0`` and ``ξ_1`` parameters dependent on the
mobile phase gas [1].

[1] Blumberg, Leonid M., Temperature-Programmed Gas Chromatography,
Wiley-VCH, 2010.
"""
function viscosity(T::Float64, gas::String)
    # using empiric model from Blumberg.2010
    if gas=="He"
        ηst = 18.63e-6
        ξ₀ = 0.6958
        ξ₁ = -0.0071
    elseif gas=="H2"
        ηst = 8.382e-6
        ξ₀ = 0.6892
        ξ₁ = 0.005
    elseif gas=="N2"
        ηst = 16.62e-6
        ξ₀ = 0.7665
        ξ₁ = -0.0378
    elseif gas=="Ar"
        ηst = 21.04e-6
        ξ₀ = 0.8131
        ξ₁ = -0.0426
    else
        error("Unknown selection of gas. Choose one of these: He, H2, N2 or Ar.")
    end
    η = ηst*(T/Tst)^(ξ₀ + ξ₁*(T-Tst)/Tst)
    return η
end

"""
    holdup_time(T, pin, pout, L, d, gas)

Calculate the hold-up time in s without a gradient.

# Arguments
* `T`: Temperature in K.
* `pin`: Inlet pressure in Pa(a).
* `pout`: Outlet pressure in Pa(g).
* `L`: Length of the capillary measured in m (meter).
* `d`: Diameter of the GC column, in m.
* `gas`: Name of the mobile phase gas.

``t_M = \\frac{128}{3}\\frac{L^2}{d^2}η\\frac{p_{in}^3-p_{out}^3}{(p_{in}^2-p_{out}^2)^2}``
"""
function holdup_time(T::Float64, pin::Float64, pout::Float64, L::Float64, d::Float64, gas::String)
    # hold-up time at temperature T (non-gradient)
	η = viscosity(T, gas)
	tM = 128/3*L^2/d^2*η*(pin^3-pout^3)/(pin^2-pout^2)^2
	return tM
end

"""
    holdup_time(t, T_itp, pin_itp, pout_itp, L, d, gas; ng=false)

Calculate the hold-up time in s at time `t` with a gradient.

# Arguments
* `t`: Time in s.
* `T_itp`: Interpolated (linear) temperature `T(x,t)`.
* `pin_itp`: Interpolated (linear) inlet pressure `pin(t)`.
* `pout_itp`: Interpolated (linear) outlet pressure `pout(t)`.
* `L`: Length of the capillary measured in m (meter).
* `d`: Diameter of the GC column, in m.
* `gas`: Name of the mobile phase gas.
* `ng`: Option to calculate the simulation without a gradient (`ng = true`,
    eq. 2)
or with a gradient (`ng = false`, eq. 1).

``t_M(t) = 64\\frac{κ_L(t)}{p_{in}(t)^2-p_{out}(t)^2} \\int_0^L
d(y)^2\\frac{p(y,t)}{T(y,t)}dy`` Eq. 1

``t_M(t) =
\\frac{128}{3}\\frac{L^2}{d^2}η(t)\\frac{p_{in}(t)^3-p_{out}(t)^3}{(p_{in}(t)^2-p_{out}(t)^2)^2}``
Eq. 2
"""
function holdup_time(t, T_itp, pin_itp, pout_itp, L, d, gas; ng=false)
    # hold-up time at time t in a temperature program with potential thermal gradient
    if ng==true
        η = viscosity(L, t, T_itp, gas)
        tM = 128/3*L^2/d(L)^2*η*(pin_itp(t)^3-pout_itp(t)^3)/(pin_itp(t)^2-pout_itp(t)^2)^2
    else
        κL = flow_restriction(L, t, T_itp, d, gas; ng=false)
        integral = quadgk(y -> d(y)^2*pressure(y, t, T_itp, pin_itp, pout_itp, L, d, gas; ng=false)/T_itp(y, t), 0, L)[1]
        tM = 64*κL/(pin_itp(t)^2-pout_itp(t)^2)*integral
    end
    return tM
end

"""
    flow(T, pin, pout, L, d, gas)

Calculate the normalized flow through the GC column in m³/s without a gradient.

# Arguments
* `T`: Temperature in K.
* `pin`: Inlet pressure in Pa(a).
* `pout`: Outlet pressure in Pa(g).
* `L`: Length of the capillary measured in m (meter).
* `d`: Diameter of the GC column, in m.
* `gas`: Name of the mobile phase gas.

``F =
\\frac{π}{256}\\frac{T_n}{p_n}\\frac{d^4}{L}\\frac{p_{in}^2-p_{out}^2}{η T}``

with ``T_n`` the normalized temperature (``T_n=(25 + 273.15)``K), ``p_n``
the normalized pressure (``p_n = 101300`` Pa(a)) and ``η`` the viscosity
the mobile phase gas at temperature ``T``.
"""
function flow(T::Float64, pin::Float64, pout::Float64, L::Float64, d::Float64, gas::String)
	# normalized Flow at temperature T (non-gradient)
	η = GasChromatographySimulator.viscosity(T, gas)
	F = π/256 * Tn/pn * d^4/L * (pin^2-pout^2)/(η*T)
	return F
end

"""
    flow(t, T_itp, pin_itp, pout_itp, L, d, gas; ng=false)

Calculate the normalized flow through the GC column in m³/s at time `t`.

# Arguments
* `t`: Time in s.
* `T_itp`: Interpolated (linear) temperature `T(x,t)`.
* `pin_itp`: Interpolated (linear) inlet pressure `pin(t)`.
* `pout_itp`: Interpolated (linear) outlet pressure `pout(t)`.
* `L`: Length of the capillary measured in m (meter).
* `d`: Diameter of the GC column, in m.
* `gas`: Name of the mobile phase gas.
* `ng`: Option to calculate the simulation without a gradient (`ng = true`,
    eq. 2)
or with a gradient (`ng = false`, eq. 1).

``F(t) =
\\frac{π}{256}\\frac{T_n}{p_n}\\frac{p_{in}(t)^2-p_{out}(t)^2}{κ_L(t)}``
Eq. 1

``F(t) =
\\frac{π}{256}\\frac{T_n}{p_n}\\frac{d^4}{L}\\frac{p_{in}(t)^2-p_{out}(t)^2}{η(t)
T(t)}``
Eq. 2

with ``T_n`` the normalized temperature (``T_n=(25 + 273.15)``K), ``p_n``
the normalized pressure (``p_n = 101300`` Pa(a)), ``κ_L`` the flow
restriction of the column and ``η`` the viscosity
the mobile phase gas at temperature ``T``.
"""
function flow(t, T_itp, pin_itp, pout_itp, L, d, gas; ng=false)
	# normalized Flow at time t in a temperature program with potential thermal
	# gradient
    # TODO: test for gradient in d(x)
    if ng==true
	    η = GasChromatographySimulator.viscosity(L, t, T_itp, gas)
	    F = π/256 * Tn/pn * d(L)^4/L * (pin_itp(t)^2-pout_itp(t)^2)/(η*T_itp(L,t))
    else
        κL = flow_restriction(L, t, T_itp, d, gas; ng=false)
        F = π/256 * Tn/pn * (pin_itp(t)^2-pout_itp(t)^2)/κL
    end
	return F
end

"""
    mobile_phase_residency(x, t, T_itp, pin_itp, pout_itp, L, d, gas; ng=false)

Calculate the residency (the inverse velocity) of the mobile phase at
position `x` at time `t`.

# Arguments
* `x`: Position along the GC column, in m.
* `t`: Time in s.
* `T_itp`: Interpolated (linear) temperature `T(x,t)`.
* `pin_itp`: Interpolated (linear) inlet pressure `pin(t)`.
* `pout_itp`: Interpolated (linear) outlet pressure `pout(t)`.
* `L`: Length of the capillary measured in m (meter).
* `d`: Diameter of the GC column, in m.
* `gas`: Name of the mobile phase gas.
* `ng`: Option to calculate the simulation without a gradient (`ng = true`)
or with a gradient (`ng = false`).

``r_M(x,t) = 64 \\frac{d^2 κ_L}{T(x,t)}\\frac{p(x,t)}{p_{in}^2-p_{out}^2}``

with ``T_n`` the normalized temperature (``T_n=(25 + 273.15)``K), ``p_n``
the normalized pressure (``p_n = 101300`` Pa(a)), ``κ_L`` the flow
restriction of the column and ``p(x,t)`` the local pressure.

See also: [`pressure`](@ref), [`flow_restriction`](@ref)
"""
function mobile_phase_residency(x, t, T_itp, pin_itp, pout_itp, L, d, gas; ng=false)
    pp = pressure(x, t, T_itp, pin_itp, pout_itp, L, d, gas; ng=ng)
    κL = flow_restriction(L, t, T_itp, d, gas; ng=ng)
    rM = 64*(pp*(d(x))^2)/T_itp(x, t)*κL/(pin_itp(t)^2-pout_itp(t)^2)
    return rM
end

"""
    residency(x, t, T_itp, pin_itp, pout_itp, L, d, df, gas, Tchar, θchar, ΔCp,  φ₀; ng=false)

Calculate the residency (the inverse velocity) of the solute at
position `x` at time `t`.

# Arguments
* `x`: Position along the GC column, in m.
* `t`: Time in s.
* `T_itp`: Interpolated (linear) temperature `T(x,t)`.
* `pin_itp`: Interpolated (linear) inlet pressure `pin(t)`.
* `pout_itp`: Interpolated (linear) outlet pressure `pout(t)`.
* `L`: Length of the capillary measured in m (meter).
* `d`: Diameter of the GC column, in m.
* `df`: Film thickness of the GC column, in m.
* `gas`: Name of the mobile phase gas.
* `Tchar`: Characteristic temperature of the solute, in K.
* `θchar`: Characteristic parameters of the solute, in °C.
* `ΔCp`: Change of the isobaric heat capacity of the solute moving from the mobile to the
stationary phase, in J mol⁻¹ K⁻¹.
* `φ₀`: Dimensionless film thickness (φ ≈ df/d) of the column for which the
thermodynamic parameters (Tchar, θchar, ΔCp) were estimated.
* `ng`: Option to calculate the simulation without a gradient (`ng = true`)
or with a gradient (`ng = false`).

``r(x,t) = r_M(x,t) \\left(1+k(x,t)\\right)``

with ``r_M`` the residency of the mobile phase and ``k(x,t)`` the retention
factor of the solute on the stationary phase.

See also: [`mobile_phase_residency`](@ref), [`retention_factor`](@ref)
"""
function residency(x, t, T_itp, pin_itp, pout_itp, L, d, df, gas, Tchar, θchar, ΔCp, φ₀; ng=false)
    r = mobile_phase_residency(x, t, T_itp, pin_itp, pout_itp, L, d, gas; ng=ng)*(1 + retention_factor(x, t, T_itp, d, df, Tchar, θchar, ΔCp, φ₀))
    return r
end

"""
    retention_factor(x, t, T_itp, d, df, Tchar, θchar, ΔCp, φ₀)

Calculate the retention factor of the solute in the stationary phase at
position `x` at time `t`.

# Arguments
* `x`: Position along the GC column, in m.
* `t`: Time in s.
* `T_itp`: Interpolated (linear) temperature `T(x,t)`.
* `d`: Diameter of the GC column, in m.
* `df`: Film thickness of the GC column, in m.
* `Tchar`: Characteristic temperature of the solute, in K.
* `θchar`: Characteristic parameters of the solute, in °C.
* `ΔCp`: Change of the isobaric heat capacity of the solute moving from the mobile to the
stationary phase, in J mol⁻¹ K⁻¹.
* `φ₀`: Dimensionless film thickness (φ ≈ df/d) of the column for which the
thermodynamic parameters (Tchar, θchar, ΔCp) were estimated.

``k(x,t) = \\frac{φ}{φ₀}
\\exp{\\left((\\frac{ΔC_p}{R}+\\frac{T_{char}}{θ_{char}})(\\frac{T_{char}}{T}+-1)
    \\frac{ΔC_p}{R}\\ln{(\\frac{T}{T_{char}})}\\right)}``

with ``R`` the molar gas constant and ``φ`` the dimensionless film thickness
of the simulated GC system (``φ = d_f/d``).

**TODO**: add option for the retention model ('ABC', 'K-centric')
"""
function retention_factor(x, t, T_itp, d, df, Tchar, θchar, ΔCp, φ₀)
    # this version of the function, where every parameter is
    # given to the function separatly seems to be the fastest
    # version
    # for now only the ideal thermodynamic model
    T = T_itp(x, t)
    φ = df(x)/d(x)
    C = ΔCp/R
    lnk₀ = (C + Tchar/θchar) * (Tchar/T - 1) + C*log(T/Tchar)
    k = φ/φ₀*exp(lnk₀)
    return k
end

"""
    plate_height(x, t, T_itp, pin_itp, pout_itp, L, d, df, gas, Tchar, θchar, ΔCp, φ₀, Dag; ng=false)

Calculate the plate height of the solute at position `x` at time `t`
according to the Golay equation.

# Arguments
* `x`: Position along the GC column, in m.
* `t`: Time in s.
* `T_itp`: Interpolated (linear) temperature `T(x,t)`.
* `pin_itp`: Interpolated (linear) inlet pressure `pin(t)`.
* `pout_itp`: Interpolated (linear) outlet pressure `pout(t)`.
* `L`: Length of the capillary measured in m (meter).
* `d`: Diameter of the GC column, in m.
* `df`: Film thickness of the GC column, in m.
* `gas`: Name of the mobile phase gas.
* `Tchar`: Characteristic temperature of the solute, in K.
* `θchar`: Characteristic parameters of the solute, in °C.
* `ΔCp`: Change of the isobaric heat capacity of the solute moving from the mobile to the
stationary phase, in J mol⁻¹ K⁻¹.
* `φ₀`: Dimensionless film thickness (φ ≈ df/d) of the column for which the
thermodynamic parameters (Tchar, θchar, ΔCp) were estimated.
* `Dag`: diffusivity of solute `a` in gas `g`.
* `ng`: Option to calculate the simulation without a gradient (`ng = true`)
or with a gradient (`ng = false`).

``H(x,t) = 2 \\frac{D_M}{u_M} + \\frac{d^2}{96}\\left(6 μ^2-16 μ +11
\\right) \\frac{u_M}{D_M} + \\frac{2}{3} d_f^2 μ(1-μ) \\frac{u_M}{D_S}``

with ``D_M`` the diffusion coefficient of the solute in the mobile phase,
``D_S`` the diffusion coefficient of the solute in the stationary phase,
``u_M`` the velocity of the mobile phase and μ the mobility of the solute.

``D_S`` is correlated to ``D_M`` by: 

``D_S = \\frac{D_M}{10000}``

**TODO**: alternative correlations?

``u_M`` is realated to the residency of the mobile phase ``r_M``:

``u_M = \\frac{1}{r_M}``

μ is correlated to the retention factor ``k``:

``μ = \\frac{1}{1 + k}``

See also: [`diffusion_mobile`](@ref), [`mobile_phase_residency`](@ref), [`retention_factor`](@ref)
"""
function plate_height(x, t, T_itp, pin_itp, pout_itp, L, d, df, gas, Tchar, θchar, ΔCp, φ₀, Dag; ng=false)
    id = d(x)# - 2.0*df(x)
    uM = 1/mobile_phase_residency(x, t, T_itp, pin_itp, pout_itp, L, d, gas; ng=ng)
    μ = 1/(1 + retention_factor(x, t, T_itp, d, df, Tchar, θchar, ΔCp, φ₀))
    DM = diffusion_mobile(x, t, T_itp, pin_itp, pout_itp, L, d, gas, Dag; ng=ng)
    DS = DM/10000
    H1 = 2*DM/uM
    H2 = id^2/96*(6*μ^2-16*μ+11)*uM/DM
    H3 = 2/3*df(x)^2*μ*(1-μ)*uM/DS
    H = H1 + H2 + H3
    return H
end

"""
    diffusion_mobile(x, t, T_itp, pin_itp, pout_itp, L, d, gas, Dag; ng=false)

Calculate the diffusion coefficient of the solute in the mobile phase at
position `x` at time `t`.

# Arguments
* `x`: Position along the GC column, in m.
* `t`: Time in s.
* `T_itp`: Interpolated (linear) temperature `T(x,t)`.
* `pin_itp`: Interpolated (linear) inlet pressure `pin(t)`.
* `pout_itp`: Interpolated (linear) outlet pressure `pout(t)`.
* `L`: Length of the capillary measured in m (meter).
* `d`: Diameter of the GC column, in m.
* `gas`: Name of the mobile phase gas.
* `Dag`: diffusivity of solute `a` in gas `g`.
* `ng`: Option to calculate the simulation without a gradient (`ng = true`)
or with a gradient (`ng = false`).

``D_M(x,t) = D_{ag} \\frac{T(x,t)^{1.75}}{p(x,t)}``
"""
function diffusion_mobile(x, t, T_itp, pin_itp, pout_itp, L, d, gas, Dag; ng=false)
    DM = T_itp(x, t)^1.75/pressure(x, t, T_itp, pin_itp, pout_itp, L, d, gas; ng=ng)*Dag
    return DM
end
#---End-Functions-of-the-physical-model---

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
		sol[i] = solving_odesystem_r(par.sys, par.prog, par.sub[i], par.opt)
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
        sol[i] = solving_migration(par.sys, par.prog, par.sub[i], par.opt)
        peak[i] = solving_peakvariance(sol[i], par.sys, par.prog, par.sub[i], par.opt)
    end
    return sol, peak
end

"""
    solving_migration(sys::System, prog::Program, sub::Substance, opt::Options)

Solve for the migration ``t(z)`` of solute `sub` in the GC system `sys` with
the program `prog` and the options `opt`.

Note: The result is the solution structure from
DifferentialEquations.jl.
"""
function solving_migration(sys::System, prog::Program, sub::Substance, opt::Options)
	f_tz(t,p,z) = residency(z, t, prog.T_itp, prog.pin_itp, prog.pout_itp, sys.L, sys.d, sys.df, sys.gas, sub.Tchar, sub.θchar, sub.ΔCp, sub.φ₀; ng=opt.ng)
    t₀ = sub.t₀
    zspan = (0.0,sys.L)
    prob_tz = ODEProblem(f_tz, t₀, zspan)
    solution_tz = solve(prob_tz, alg=opt.alg, abstol=opt.abstol,reltol=opt.reltol)
    return solution_tz
end

"""
    solving_peakvariance(solution_tz, sys::System, prog::Program, sub::Substance, opt::Options)

Solve for the development of the peak variance ``τ²(z)`` of solute `sub` in the GC system `sys` with
the program `prog` and the options `opt` during its migration defined by `solution_tz`.

Note: The result is the solution structure from
DifferentialEquations.jl.
"""
function solving_peakvariance(solution_tz, sys, prog, sub, opt)
    t(z) = solution_tz(z)
    p = [sys, prog, sub, opt]
    f_τ²z(τ²,p,z) = peakode(z, t(z), τ², sys, prog, sub, opt)
    τ²₀ = sub.τ₀^2
    zspan = (0.0, sys.L)
    prob_τ²z = ODEProblem(f_τ²z, τ²₀, zspan, p)
    solution_τ²z = solve(prob_τ²z, alg=opt.alg, abstol=opt.abstol,reltol=opt.reltol)
    return solution_τ²z
end

"""
    solving_odesystem_r(sys::System, prog::Program, sub::Substance, opt::Options)

Solve the migration ``t(z)`` and peak variance development ``τ²(z)`` of solute `sub` in the GC system `sys` with
the program `prog` and the options `opt` as a system of ODEs.

Note: The result is the solution structure from
DifferentialEquations.jl.

See also: [`odesystem_r!`](@ref)
"""
function solving_odesystem_r(sys::System, prog::Program, sub::Substance, opt::Options)
    t₀ = [sub.t₀; sub.τ₀^2]
    zspan = (0.0,sys.L)
	p = [sys, prog, sub, opt]
    prob = ODEProblem(odesystem_r!, t₀, zspan, p)

    solution = solve(prob, alg=opt.alg, abstol=opt.abstol,reltol=opt.reltol)

    if solution.t[end]<sys.L
        solution = solve(prob, alg=opt.alg, abstol=opt.abstol,reltol=opt.reltol, dt=sys.L/1000000)
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
	sys = p[1]
	prog = p[2]
	sub = p[3]
	opt = p[4]
    dt[1] = residency(z, t[1], prog.T_itp, prog.pin_itp, prog.pout_itp, sys.L, sys.d, sys.df, sys.gas, sub.Tchar, sub.θchar, sub.ΔCp, sub.φ₀; ng=opt.ng)
    dt[2] = peakode(z, t[1], t[2], sys, prog, sub, opt)
end

"""
    peakode(z, t, τ², sys, prog, sub, opt)

The second ODE function for the ODE system describing the peak variance
development ``τ²(z)``, using (in parts) automatic differentiation.

``\\frac{dτ^2}{dz} = H(z, t(z)) r(z, t(z)) + 2 τ^2(z, t(z))
\\frac{∂r}{∂t}(z,t(z))``

**TODO**: alternative to QuadGK.jl for integration which is available for
ForwardDiff.jl 

See also: [`solving_odesystem_r`](@ref), [`odesystem_r!`](@ref)
"""
function peakode(z, t, τ², sys, prog, sub, opt)
	# alternative function
    if opt.ng==true
        r_ng(zt) = residency(zt[1], zt[2], prog.T_itp, prog.pin_itp, prog.pout_itp, sys.L, sys.d, sys.df, sys.gas, sub.Tchar, sub.θchar, sub.ΔCp, sub.φ₀; ng=true)
        H_ng(z,t) = plate_height(z, t, prog.T_itp, prog.pin_itp, prog.pout_itp, sys.L, sys.d, sys.df, sys.gas, sub.Tchar, sub.θchar, sub.ΔCp, sub.φ₀, sub.Dag; ng=true)
        ∂r∂t_ng(z,t) = ForwardDiff.gradient(r_ng, [z, t])[2]
        return H_ng(z,t)*r_ng([z,t])^2 + 2*τ²*∂r∂t_ng(z,t)
    else
		d(z) = sys.d(z)
		T(z,t) = prog.T_itp(z, t)
		η(z,t) = viscosity(z, t, prog.T_itp, sys.gas)
		c(zt) = η(zt[1], zt[2])*T(zt[1], zt[2])
		∂c∂t(z,t) = ForwardDiff.gradient(c, [z,t])[2]
		pi2(t) = prog.pin_itp(t)^2
		po2(t) = prog.pout_itp(t)^2
		∂pi2∂t(t) = ForwardDiff.derivative(pi2, t)
		∂po2∂t(t) = ForwardDiff.derivative(po2, t)
        κ(z,t) = flow_restriction(z, t, prog.T_itp, sys.d, sys.gas)
		κL(t) = κ(sys.L,t)
		∂κ∂t(z,t) = quadgk(y -> d(y)^-4*∂c∂t(y,t), 0.0, z)[1]
        ∂κL∂t(t) = ∂κ∂t(sys.L,t)
		e(z,t) = κ(z,t)/κL(t)
		∂e∂t(z,t) = (∂κ∂t(z,t)*κL(t)-κ(z,t)*∂κL∂t(t))/κL(t)^2
		p(z,t) = pressure(z, t, prog.T_itp, prog.pin_itp, prog.pout_itp, sys.L, sys.d, sys.gas)
		∂p∂t(z,t) = 1/(2*p(z,t))*(∂pi2∂t(t)-(∂e∂t(z,t)*(pi2(t)-po2(t))+e(z,t)*(∂pi2∂t(t)-∂po2∂t(t))))
		rM(z,t) = mobile_phase_residency(z, t, prog.T_itp, prog.pin_itp, prog.pout_itp, sys.L, sys.d, sys.gas)
		a(z,t) = κL(t)*p(z,t)
		∂a∂t(z,t) = ∂κL∂t(t)*p(z,t)+κL(t)*∂p∂t(z,t)
		b(zt) = T(zt[1],zt[2])*(pi2(zt[2])-po2(zt[2]))
		∂b∂t(z,t) = ForwardDiff.gradient(b, [z,t])[2]
		∂rM∂t(z,t) = 64*d(z)^2*((∂a∂t(z,t)*b([z,t])-a(z,t)*∂b∂t(z,t))/b([z,t])^2)
		
		k(zt) = retention_factor(zt[1], zt[2], prog.T_itp, sys.d, sys.df, sub.Tchar, sub.θchar, sub.ΔCp, sub.φ₀)
        ∂k∂t(z,t) = ForwardDiff.gradient(k, [z, t])[2]
		
		r(z,t) = residency(z, t, prog.T_itp, prog.pin_itp, prog.pout_itp, sys.L, sys.d, sys.df, sys.gas, sub.Tchar, sub.θchar, sub.ΔCp, sub.φ₀)
        H(z,t) = plate_height(z, t, prog.T_itp, prog.pin_itp, prog.pout_itp, sys.L, sys.d, sys.df, sys.gas, sub.Tchar, sub.θchar, sub.ΔCp, sub.φ₀, sub.Dag)
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
        if sol[i].t[end]==par.sys.L
            tR[i] = sol[i].u[end][1]
            TR[i] = par.prog.T_itp(par.sys.L, tR[i]) - 273.15 
            uR[i] = 1/residency(par.sys.L, tR[i], par.prog.T_itp, par.prog.pin_itp, par.prog.pout_itp, par.sys.L, par.sys.d, par.sys.df, par.sys.gas, par.sub[i].Tchar, par.sub[i].θchar, par.sub[i].ΔCp, par.sub[i].φ₀)
            τR[i] = sqrt(sol[i].u[end][2])
            σR[i] = τR[i]*uR[i]
            kR[i] = retention_factor(par.sys.L, tR[i], par.prog.T_itp, par.sys.d, par.sys.df, par.sub[i].Tchar, par.sub[i].θchar, par.sub[i].ΔCp, par.sub[i].φ₀)
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
        if sol[i].t[end]==par.sys.L
            tR[i] = sol[i].u[end]
            TR[i] = par.prog.T_itp(par.sys.L, tR[i]) - 273.15 
            uR[i] = 1/residency(par.sys.L, tR[i], par.prog.T_itp, par.prog.pin_itp, par.prog.pout_itp, par.sys.L, par.sys.d, par.sys.df, par.sys.gas, par.sub[i].Tchar, par.sub[i].θchar, par.sub[i].ΔCp, par.sub[i].φ₀)
            τR[i] = sqrt(peak[i].u[end])
            σR[i] = τR[i]*uR[i]
            kR[i] = retention_factor(par.sys.L, tR[i], par.prog.T_itp, par.sys.d, par.sys.df, par.sub[i].Tchar, par.sub[i].θchar, par.sub[i].ΔCp, par.sub[i].φ₀)
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
#---End-Result-Functions---

end # module
