module GasChromatographySimulator

using Reexport
using Interpolations
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

# structures
"""
    System(L, d, a_d, df, a_df, sp, gas)

    Structure describing the GC system. The system consists of a
    capillary (called column) of length `L` and inner diameter `d`. The capillary is coated
    with thin film (thickness `df`) of a stationary phase `sp`. A mobile phase
    consisting of a `gas` moves through the capillary. The diameter `d` and film
    thickness `d_f` are described as functions to enable the option of gradients of
    diameter or film thickness.

    # Arguments
    * `L`: Length of the capillary measured in m (meter)
    * `d`: A function `d(x, a_d)` of `x`, the position along the capillary,
    describing the diameter in m (meter).
    * `a_d`: Parameters of the diameter function. 
    * `d_f`: A function `d_f(x, a_df)` of `x`, describing the film thickness in m
    (meter).
    * `sp`: The name of the stationary phase.
    * `gas`: THe name of the mobile phase. Allowed values: He, H2 or N2.
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

    Description of the setup ...

    # Arguments
    * `time_steps`: Time steps in s (seconds). 
    * `temp_steps`: Temperature steps in °C (degree celsius).
    * `pin_steps`: Inlet pressure steps in Pa(a) (pascal, absolute).
    * `pout_steps`: Outlet pressure steps in Pa(a) (pascal, absolute).
    * `gf`: Gradient function `gf(x, a_gf)`, describes the thermal gradient.
    * `a_gf`: Parameters of the gradient function.
    * `T_itp`: Interpolated (linear) temperature `T(x,t)`, constructed from
      `time_steps`, `temp_steps` and `gf`.
    * `pin_itp`: Interpolated (linear) inlet pressure `pin(t)`, constructed from
      `time_steps` and `pin_steps`.
    * `pout_itp`: Interpolated (linear) outlet pressure `pout(t)`, constructed
      from `time_steps` and `pout_steps`.  
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

Base.@kwdef struct Options
    alg                 # algorithmen for the ODE solver
    abstol              # absolute tolerance for ODE solver
    reltol              # relative tolerance for ODE solver 
    Tcontrol::String    # temperature control at 'inlet' (top) or 'outlet' (bottom) of the column
	odesys::Bool  		# calculate the two ODEs (migration and peak-width) separately (false) or 
                        # combined as a system of ODEs (true)
end

# test
struct Parameters
    sys::System
    prog::Program
    sub::Array{Substance,1}
    opt::Options
end

function constructor_System(L, d, df, sp, gas)
    # function to construct the System structure
    # for the case of constant diameter and constant film thickness
    d_func(x) = gf_const(x, [d])
    df_func(x) = gf_const(x, [df])
    sys = System(L, d_func, [d], df_func, [df], sp, gas)
    return sys
end

function constructor_Program(time_steps::Array{<:Real, 1}, temp_steps::Array{<:Real, 1}, pin_steps::Array{<:Real, 1}, pout_steps::Array{<:Real, 1}, ΔT_steps::Array{<:Real, 1}, x₀_steps::Array{<:Real, 1}, L₀_steps::Array{<:Real, 1}, α_steps::Array{<:Real, 1}, Tcontrol, L)
    # function to construct the Program structure
    # using as gradient function the exponential model 'gf_exp(x,a_gf,Tcontrol)'
    a_gf = [ΔT_steps x₀_steps L₀_steps α_steps]
    gf(x) = gf_exp(x, a_gf, Tcontrol)
    T_itp = temperature_interpolation(time_steps, temp_steps, gf, L)
    pin_itp = pressure_interpolation(time_steps, pin_steps)
    pout_itp = pressure_interpolation(time_steps, pout_steps)
    prog = Program(time_steps, temp_steps, pin_steps, pout_steps, gf, a_gf, T_itp, pin_itp, pout_itp)
    return prog
end

function constructor_Program(time_steps::Array{<:Real, 1}, temp_steps::Array{<:Real, 1}, pin_steps::Array{<:Real, 1}, pout_steps::Array{<:Real, 1}, L)
    # function to construct the Program structure
    # without a thermal gradient
    # using as gradient function the exponential model 'gf_exp(x,a_gf,Tcontrol)'
    a_gf = [zeros(length(time_steps)) zeros(length(time_steps)) L.*ones(length(time_steps)) zeros(length(time_steps))]
    gf(x) = gf_exp(x, a_gf, "inlet")
    T_itp = temperature_interpolation(time_steps, temp_steps, gf, L)
    pin_itp = pressure_interpolation(time_steps, pin_steps)
    pout_itp = pressure_interpolation(time_steps, pout_steps)
    prog = Program(time_steps, temp_steps, pin_steps, pout_steps, gf, a_gf, T_itp, pin_itp, pout_itp)
    return prog
end

# gradient functions
function gf_const(x,a::Array{<:Real,1})
    # fixed parameter values
    f = a[1]
    return f
end

function gf_const(x,a::Array{<:Real,2})
    # parameters can change (in time)
    f = a[:,1]
    return f
end

function gf_linear(x,a::Array{<:Real,1})
    # fixed parameter values
    f₀ = a[1] # start value
    x₀ = a[2] # shift in x 
    L₀ = a[3] # length of the linear segment
    # gradient is -f₀/L₀
    f = f₀ * (1 - (x-x₀) / L₀)
end

function gf_linear(x,a::Array{<:Real,2})
    # parameters can change (in time)
    f₀ = a[:,1] # start value
    x₀ = a[:,2] # shift in x 
    L₀ = a[:,3] # length of the linear segment
    # gradient is -f₀/L₀
    f = f₀ .* (1 .- (x.-x₀) ./ L₀)
end

function gf_exp(x,a::Array{<:Real,1})
    # fixed parameter values
    f₀ = a[1] # start value
    x₀ = a[2] # shift in x 
    L₀ = a[3] # length of the linear segment
    α  = a[4] # exponential factor, α=0 -> linear
    f = f₀ .* (1 .- exp.(α.*(1 .- (x.-x₀) ./ L₀)) + (1 .- (x.-x₀) ./ L₀) .* exp.(α))
end

function gf_exp(x,a::Array{<:Real,2}, Tcontrol::String)
    # parameters can change (in time)
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
    return f
end

function temperature_interpolation(time_steps::Array{<:Real,1}, temp_steps::Array{<:Real,1}, gradient_function::Function, L::Float64)
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

function pressure_interpolation(time_steps::Array{<:Real,1}, press_steps::Array{<:Real,1})
    p_itp = LinearInterpolation((cumsum(time_steps), ), press_steps, extrapolation_bc=Flat())
    return p_itp
end

# test
function load_solute_database(db_path::String, db::String, sp::String, gas::String, solutes::Array{String,1}, t₀::Array{Float64,1}, τ₀::Array{Float64,1})
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
	#
	# at the moment only the new database config is supported (new data is appended in new rows)
	# also, if for a solute no data is available in the database no error is given!!!
	df_solute_db = DataFrame(CSV.File(string(db_path,"/",db), header=1, silencewarnings=true))
	if size(df_solute_db)[2]>14
		# old database config (additional stationary phases are appended in new columns, only one row for a solute)
        error("Data format not supported. Use the appended database structure.")
	else
		# new database config (simpler, additional stationary phases are appended in new rows, multiple entrys for a solute)
		# 1. Filter the stationary phase
		phases_db = unique(df_solute_db.Phase)
		if isa(findfirst(phases_db.==sp), Nothing) && sp!=""
			error("Unknown selction of stationary phase. Choose one of these: $phases_db")
		elseif sp=="" # no stationary phase is selected, like for transferlines
			df_solute_db_filtered = df_solute_db[!,1:8] # use only the data of the first 8 columns
			# 2. Filter the solutes, not using multiple entrys
			df_solute_db_filtered_1=unique(df_solute_db_filtered[in(solutes).(df_solute_db_filtered.Name),:])
			# use placeholder values
			Tchar = ones(length(solutes))
			θchar = ones(length(solutes))
			ΔCp = ones(length(solutes))
			φ₀ = ones(length(solutes))
			Annotation = fill("no sp", length(solutes))
		else
			df_solute_db_filtered = filter([:Phase] => x -> x==sp, df_solute_db)
			# 2. Filter the solutes
			df_solute_db_filtered_1=df_solute_db_filtered[in(solutes).(df_solute_db_filtered.Name),:]
			# values
			Tchar = df_solute_db_filtered_1.Tchar.+Tst
			θchar = df_solute_db_filtered_1.thetachar
			ΔCp = df_solute_db_filtered_1.DeltaCp
			φ₀ = df_solute_db_filtered_1.phi0
			Annotation = df_solute_db_filtered_1.Annotation
		end
		# 3. Construct the sub()-structure
		sub = Array{Substance}(undef, size(df_solute_db_filtered_1)[1])
		for i=1:size(df_solute_db_filtered_1)[1]
			Dag = diffusivity(df_solute_db_filtered_1.Molmass[i], 
											df_solute_db_filtered_1.Cnumber[i], 
											df_solute_db_filtered_1.Hnumber[i], 
											df_solute_db_filtered_1.Onumber[i], 
											df_solute_db_filtered_1.Nnumber[i], 
											df_solute_db_filtered_1.Ringnumber[i], 
											gas)
			sub[i] = Substance(df_solute_db_filtered_1.Name[i],
										df_solute_db_filtered_1.CAS[i],
										Tchar[i], 
										θchar[i], 
										ΔCp[i], 
										φ₀[i],
										Annotation[i],
										Dag, 
										τ₀[i],
										t₀[i])
		end
	end
	return sub
end

function diffusivity(M::Float64, Cn::Real, Hn::Real, On::Real, Nn::Real, Rn::Real, gas::String)
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

# functions of the physical model
function pressure(x, t, T_itp, pin_itp, pout_itp, L, d, gas; ng=false)
    if ng==true
        pp = sqrt(pin_itp(t)^2 - x/L*(pin_itp(t)^2-pout_itp(t)^2))
    else
        pp = sqrt(pin_itp(t)^2 - flow_restriction(x, t, T_itp, d, gas)/flow_restriction(L, t, T_itp, d, gas)*(pin_itp(t)^2-pout_itp(t)^2))
    end
    return pp
end

function flow_restriction(x, t, T_itp, d, gas; ng=false)
    if ng==true
        κ = x*viscosity(x, t, T_itp, gas)*T_itp(x, t)*d(x)^-4
    else
        κ = quadgk(y -> viscosity(y, t, T_itp, gas)*T_itp(y, t)*d(y)^-4, 0, x)[1]
    end
    return κ
end

function viscosity(x, t, T_itp, gas::String)
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

function mobile_phase_residency(x, t, T_itp, pin_itp, pout_itp, L, d, gas; ng=false)
    pp = pressure(x, t, T_itp, pin_itp, pout_itp, L, d, gas; ng=ng)
    κL = flow_restriction(L, t, T_itp, d, gas; ng=ng)
    rM = 64*(pp*(d(x))^2)/T_itp(x, t)*κL/(pin_itp(t)^2-pout_itp(t)^2)
    return rM
end

function residency(x, t, T_itp, pin_itp, pout_itp, L, d, df, gas, ΔCp, Tchar, θchar, φ₀; ng=false)
    r = mobile_phase_residency(x, t, T_itp, pin_itp, pout_itp, L, d, gas; ng=ng)*(1 + retention_factor(x, t, T_itp, d, df, ΔCp, Tchar, θchar, φ₀))
    return r
end

function retention_factor(x, t, T_itp, d, df, ΔCp, Tchar, θchar, φ₀)
    # this version of the function, where every parameter is
    # given to the function separatly seems to be the fastest
    # version
    # for now only the ideal thermodynamic model
    T = T_itp(x, t)
    φ::Float64 = df(x)/d(x)
    C = ΔCp/R
    lnk₀ = (C + Tchar/θchar) * (Tchar/T - 1) + C*log(T/Tchar)
    k = φ/φ₀*exp(lnk₀)
    return k
end

function plate_height(x, t, T_itp, pin_itp, pout_itp, L, d, df, gas, ΔCp, Tchar, θchar, φ₀, Dag; ng=false)
    id = d(x)# - 2.0*df(x)
    uM = 1/mobile_phase_residency(x, t, T_itp, pin_itp, pout_itp, L, d, gas; ng=ng)
    μ = 1/(1 + retention_factor(x, t, T_itp, d, df, ΔCp, Tchar, θchar, φ₀))
    DM = diffusion_mobile(x, t, T_itp, pin_itp, pout_itp, L, d, gas, Dag; ng=ng)
    DS = DM/10000
    H1 = 2*DM/uM
    H2 = id^2/96*(6*μ^2-16*μ+11)*uM/DM
    H3 = 2/3*df(x)^2*μ*(1-μ)*uM/DS
    H = H1 + H2 + H3
    return H
end

function diffusion_mobile(x, t, T_itp, pin_itp, pout_itp, L, d, gas, Dag; ng=false)
    DM = T_itp(x, t)^1.75/pressure(x, t, T_itp, pin_itp, pout_itp, L, d, gas; ng=ng)*Dag
    return DM
end

# solving functions
# test
function solve_system_multithreads(par; ng=false)
	n = length(par.sub)
	sol = Array{Any}(undef, n)
	Threads.@threads for i=1:n
		sol[i] = solving_odesystem_r(par.sys, par.prog, par.sub[i], par.opt; ng=ng)
	end
	return sol
end

# test
function solve_multithreads(par; ng=false)
    n = length(par.sub)
    sol = Array{Any}(undef, n)
    peak = Array{Any}(undef, n)
    Threads.@threads for i=1:n
        sol[i] = solving_migration(par.sys, par.prog, par.sub[i], par.opt; ng=ng)
        peak[i] = solving_peakvariance(sol[i], par.sys, par.prog, par.sub[i], par.opt; ng=ng)
    end
    return sol, peak
end

function solving_migration(sys, prog, sub, opt; ng=false)
    #---------------------------------------------------------------------------
    # solves the differential equations for migration t(z) of a solute with the
    # parameters p
    #---------------------------------------------------------------------------
	f_tz(t,p,z) = residency(z, t, prog.T_itp, prog.pin_itp, prog.pout_itp, sys.L, sys.d, sys.df, sys.gas, sub.ΔCp, sub.Tchar, sub.θchar, sub.φ₀; ng=ng)
    t₀ = sub.t₀
    zspan = (0.0,sys.L)
    prob_tz = ODEProblem(f_tz, t₀, zspan)
    solution_tz = solve(prob_tz, alg=opt.alg, abstol=opt.abstol,reltol=opt.reltol)
    return solution_tz
end

function solving_peakvariance(solution_tz, sys, prog, sub, opt; ng=false)
    #---------------------------------------------------------------------------
    # solves the differential equations for development of the peak variance
    # σ²(z) of a solute with migration solution solution_tz and the
    # parameters p
    # using automatic differentiation ForwardDiff and manuall differentiation
    #---------------------------------------------------------------------------
    t(z) = solution_tz(z)
    p = [sys, prog, sub, opt]
    f_τ²z(τ²,p,z) = peakode(z, t(z), τ², sys, prog, sub, opt; ng=ng)
    τ²₀ = sub.τ₀^2
    zspan = (0.0, sys.L)
    prob_τ²z = ODEProblem(f_τ²z, τ²₀, zspan, p)
    solution_τ²z = solve(prob_τ²z, alg=opt.alg, abstol=opt.abstol,reltol=opt.reltol)
    return solution_τ²z
end

function solving_odesystem_r(sys, prog, sub, opt; ng=false)
    #---------------------------------------------------------------------------
    # solves the differential equations for migration t(z) of a solute with the
    # parameters p
    #---------------------------------------------------------------------------
    t₀ = [sub.t₀; sub.τ₀^2]
    zspan = (0.0,sys.L)
	p = [sys, prog, sub, opt, ng]
    prob = ODEProblem(odesystem_r!, t₀, zspan, p)
    #try
        solution = solve(prob, alg=opt.alg, abstol=opt.abstol,reltol=opt.reltol)
    #catch
    #    solution = solve(prob, alg=p.opt.alg, abstol=p.opt.abstol,reltol=p.opt.reltol, dt=p.sys.L/1000)
    #end
    if solution.t[end]<sys.L
        solution = solve(prob, alg=opt.alg, abstol=opt.abstol,reltol=opt.reltol, dt=sys.L/1000000)
    end
    return solution
end

function odesystem_r!(dt, t, p, z)
    # t[1] ... t time
    # t[2] ... τ² band variance
	sys = p[1]
	prog = p[2]
	sub = p[3]
	opt = p[4]
    ng = p[5]
    dt[1] = residency(z, t[1], prog.T_itp, prog.pin_itp, prog.pout_itp, sys.L, sys.d, sys.df, sys.gas, sub.ΔCp, sub.Tchar, sub.θchar, sub.φ₀; ng=ng)
    dt[2] = peakode(z, t[1], t[2], sys, prog, sub, opt; ng=ng)
end

function peakode(z, t, τ², sys, prog, sub, opt; ng=false)
	# alternative function
    if ng==true
        r_ng(zt) = residency(zt[1], zt[2], prog.T_itp, prog.pin_itp, prog.pout_itp, sys.L, sys.d, sys.df, sys.gas, sub.ΔCp, sub.Tchar, sub.θchar, sub.φ₀; ng=true)
        H_ng(z,t) = plate_height(z, t, prog.T_itp, prog.pin_itp, prog.pout_itp, sys.L, sys.d, sys.df, sys.gas, sub.ΔCp, sub.Tchar, sub.θchar, sub.φ₀, sub.Dag; ng=true)
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
		
		k(zt) = retention_factor(zt[1], zt[2], prog.T_itp, sys.d, sys.df, sub.ΔCp, sub.Tchar, sub.θchar, sub.φ₀)
        ∂k∂t(z,t) = ForwardDiff.gradient(k, [z, t])[2]
		
		r(z,t) = residency(z, t, prog.T_itp, prog.pin_itp, prog.pout_itp, sys.L, sys.d, sys.df, sys.gas, sub.ΔCp, sub.Tchar, sub.θchar, sub.φ₀)
        H(z,t) = plate_height(z, t, prog.T_itp, prog.pin_itp, prog.pout_itp, sys.L, sys.d, sys.df, sys.gas, sub.ΔCp, sub.Tchar, sub.θchar, sub.φ₀, sub.Dag)
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

# examine the results
# test
function peaklist(sol, p)
	n = length(p.sub)
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
        Name[i] = p.sub[i].name
        if sol[i].t[end]==p.sys.L
            tR[i] = sol[i].u[end][1]
            TR[i] = p.prog.T_itp(p.sys.L, tR[i]) - 273.15 
            uR[i] = 1/residency(p.sys.L, tR[i], p.prog.T_itp, p.prog.pin_itp, p.prog.pout_itp, p.sys.L, p.sys.d, p.sys.df, p.sys.gas, p.sub[i].ΔCp, p.sub[i].Tchar, p.sub[i].θchar, p.sub[i].φ₀)
            τR[i] = sqrt(sol[i].u[end][2])
            σR[i] = τR[i]*uR[i]
            kR[i] = retention_factor(p.sys.L, tR[i], p.prog.T_itp, p.sys.d, p.sys.df, p.sub[i].ΔCp, p.sub[i].Tchar, p.sub[i].θchar, p.sub[i].φ₀)
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

# test
function peaklist(sol, peak, p)
	n = length(p.sub)
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
        Name[i] = p.sub[i].name
        if sol[i].t[end]==p.sys.L
            tR[i] = sol[i].u[end]
            TR[i] = p.prog.T_itp(p.sys.L, tR[i]) - 273.15 
            uR[i] = 1/residency(p.sys.L, tR[i], p.prog.T_itp, p.prog.pin_itp, p.prog.pout_itp, p.sys.L, p.sys.d, p.sys.df, p.sys.gas, p.sub[i].ΔCp, p.sub[i].Tchar, p.sub[i].θchar, p.sub[i].φ₀)
            τR[i] = sqrt(peak[i].u[end])
            σR[i] = τR[i]*uR[i]
            kR[i] = retention_factor(p.sys.L, tR[i], p.prog.T_itp, p.sys.d, p.sys.df, p.sub[i].ΔCp, p.sub[i].Tchar, p.sub[i].θchar, p.sub[i].φ₀)
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

# test
function sol_extraction(sol, p)
    # extract the points z=t, t=u1, τ²=u2 from the solution of
    # the ODE system
	n = length(p.sub)
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
        solutes[i] = p.sub[i].name
    end
    df_sol = DataFrame(name=solutes, z=sol_z, t=sol_t, τ²=sol_τ²)
    return df_sol
end

# test
function sol_extraction(sol_tz, peak_τz, p)
    # extract the points z=t, t=u, from the solution of
    # the first ODE (sol_tz)
    # and the points z=t, τ²=u, from the solution of 
    # the second ODE (peak_τz)
	n = length(p.sub)
    sol_z = Array{Any}(undef, n)
    sol_t = Array{Any}(undef, n)
    peak_z = Array{Any}(undef, n)
    peak_τ² = Array{Any}(undef, n)
    solutes = Array{String}(undef, n)
    for i=1:n
        sol_z[i] = sol_tz[i].t
        sol_t[i] = sol_tz[i].u
        peak_z[i] = peak_τz[i].t
        peak_τ²[i] = peak_τz[i].u
		solutes[i] = p.sub[i].name
	end
    df_sol = DataFrame(name=solutes, z_t=sol_z, t=sol_t, z_τ²=peak_z, τ²=peak_τ²)
    return df_sol
end

end # module
