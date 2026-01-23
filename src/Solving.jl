#---Begin-Solving-Functions---
"""
    simulate(par::Parameters; kwargs...)

Simulate the GC system defined by the structure `par`.

Alternative call:
`simulate(L, d, df, gas, T_itp, Fpin_itp, pout_itp,Name, CAS, ann, Tchar, θchar, ΔCp, φ₀, Cag, t₀, τ₀, opt; kwargs...)`

Note: Based on the option for `odesys` the result is different. For `odesys =
true` the result is a dataframe (the peaklist) and the solution of the ODEs
as a system (solution structure from DifferentialEquations.jl). If `odesys =
false` the result is a dataframe (the peaklist) and the two solutions of the
ODEs for time ``t(z)`` and peak variance ``τ²(z)``. `kwargs...` to pass additional options to the ODE solve function as named tuples.
"""
function simulate(par; kwargs...)
    if par.opt.odesys==true
        sol = solve_system_multithreads(par; kwargs...)
    	pl = GasChromatographySimulator.peaklist(sol, par)
        return pl, sol
	else
		sol, peak = solve_separate_multithreads(par; kwargs...)
    	pl = GasChromatographySimulator.peaklist(sol, peak, par)
        return pl, sol, peak
	end
end

function simulate(L, d, df, gas, T_itp, Fpin_itp, pout_itp, Name, CAS, ann, Tchar, θchar, ΔCp, φ₀, Cag, t₀, τ₀, opt; kwargs...)
    if opt.odesys==true
        sol = solve_system_multithreads(L, d, df, gas, T_itp, Fpin_itp, pout_itp, Tchar, θchar, ΔCp, φ₀, Cag, t₀, τ₀, opt; kwargs...)
    	# TODO: either create par from the input variables 
        # or create new peaklist function
        # perhaps add some parameters like 'Name' of the solutes
        pl = GasChromatographySimulator.peaklist(sol, L, d, df, gas, T_itp, Fpin_itp, pout_itp, Name, CAS, ann, Tchar, θchar, ΔCp, φ₀, opt)
        return pl, sol
	else
		sol, peak = solve_separate_multithreads(L, d, df, gas, T_itp, Fpin_itp, pout_itp, Tchar, θchar, ΔCp, φ₀, Cag, t₀, τ₀, opt; kwargs...)
    	pl = GasChromatographySimulator.peaklist(sol, peak, L, d, df, gas, T_itp, Fpin_itp, pout_itp, Name, CAS, ann, Tchar, θchar, ΔCp, φ₀, opt)
        return pl, sol, peak
	end
end

"""
    solve_system_multithreads(par::Parameters, kwargs...)

Simulate the GC system defined by the structure `par` by solving the
ODEs for ``t(z)`` and ``τ²(z)`` together as a system of ODEs using multiple
threads (parallel computing) for the simulation of different solutes. `kwargs...` to pass additional options to the ODE solve function as named tuples. 

Note: The result is an array of the solution structure from DifferentialEquations.jl.

Alternative call:
`solve_system_multithreads(L, d, df, gas, T_itp, Fpin_itp, pout_itp, Tchar, θchar, ΔCp, φ₀, Cag, t₀, τ₀, opt; kwargs...)`
with the substance realted quantitites beeing vectors.

# Examples

```julia
julia> sol = solve_system_multithreads(par)
```
"""
function solve_system_multithreads(par; kwargs...)
	n = length(par.sub)
	sol = Array{Any}(undef, n)
	Threads.@threads for i=1:n
		sol[i] = solving_odesystem_r(par.col, par.prog, par.sub[i], par.opt; kwargs...)
	end
	return sol
end

function solve_system_multithreads(L, d, df, gas, T_itp, Fpin_itp, pout_itp, Tchar, θchar, ΔCp, φ₀, Cag, t₀, τ₀, opt; kwargs...)
	n = length(Tchar)
	sol = Array{Any}(undef, n)
	Threads.@threads for i=1:n
		sol[i] = solving_odesystem_r(L, d, df, gas, T_itp, Fpin_itp, pout_itp, Tchar[i], θchar[i], ΔCp[i], φ₀[i], Cag[i], t₀[i], τ₀[i], opt; kwargs...)
	end
	return sol
end

"""
    solve_separate_multithreads(par::Parameters, kwargs...)

Simulate the GC system defined by the structure `par` by solving the
ODEs for ``t(z)`` and ``τ²(z)`` separatly (solving ``t(z)`` and using this result
to solve for ``τ²(z)``) using multiple threads (parallel computing) for the
simulation of different solutes. `kwargs...` to pass additional options to the ODE solve function as named tuples.

Note: The result are two arrays of the solution structure from
DifferentialEquations.jl.

Alternative call:
`solve_separate_multithreads(L, d, df, gas, T_itp, Fpin_itp, pout_itp, Tchar, θchar, ΔCp, φ₀, Cag, t₀, τ₀, opt; kwargs...)`
with the substance realted quantitites beeing vectors.

# Examples

```julia
julia> sol, peak = solve_separate_multithreads(par)
```
"""
function solve_separate_multithreads(par; kwargs...)
    n = length(par.sub)
    sol = Array{Any}(undef, n)
    peak = Array{Any}(undef, n)
    Threads.@threads for i=1:n
        sol[i] = solving_migration(par.col, par.prog, par.sub[i], par.opt; kwargs...)
        peak[i] = solving_peakvariance(sol[i], par.col, par.prog, par.sub[i], par.opt; kwargs...)
    end
    return sol, peak
end

function solve_separate_multithreads(L, d, df, gas, T_itp, Fpin_itp, pout_itp, Tchar, θchar, ΔCp, φ₀, Cag, t₀, τ₀, opt; kwargs...)
    n = length(Tchar)
    sol = Array{Any}(undef, n)
    peak = Array{Any}(undef, n)
    Threads.@threads for i=1:n
        sol[i] = solving_migration(L, d, df, gas, T_itp, Fpin_itp, pout_itp, Tchar[i], θchar[i], ΔCp[i], φ₀[i], t₀[i], opt; kwargs...)
        peak[i] = solving_peakvariance(sol[i], L, d, df, gas, T_itp, Fpin_itp, pout_itp, Tchar[i], θchar[i], ΔCp[i], φ₀[i], Cag[i], τ₀[i], opt; kwargs...)
    end
    return sol, peak
end

"""
    solve_system(par::Parameters, kwargs...)

Simulate the GC system defined by the structure `par` by solving the
ODEs for ``t(z)`` and ``τ²(z)`` together as a system of ODEs using multiple
threads (parallel computing) for the simulation of different solutes. `kwargs...` to pass additional options to the ODE solve function as named tuples. . No multi-threads are used.

Note: The result is an array of the solution structure from DifferentialEquations.jl.

Alternative call:
`solve_system(L, d, df, gas, T_itp, Fpin_itp, pout_itp, Tchar, θchar, ΔCp, φ₀, Cag, t₀, τ₀, opt; kwargs...)`
with the substance realted quantitites beeing vectors.

# Examples

```julia
julia> sol = solve_system(par)
```
"""
function solve_system(par; kwargs...)
	n = length(par.sub)
	sol = Array{Any}(undef, n)
	for i=1:n
		sol[i] = solving_odesystem_r(par.col, par.prog, par.sub[i], par.opt; kwargs...)
	end
	return sol
end

function solve_system(L, d, df, gas, T_itp, Fpin_itp, pout_itp, Tchar, θchar, ΔCp, φ₀, Cag, t₀, τ₀, opt; kwargs...)
	n = length(Tchar)
	sol = Array{Any}(undef, n)
	for i=1:n
		sol[i] = solving_odesystem_r(L, d, df, gas, T_itp, Fpin_itp, pout_itp, Tchar[i], θchar[i], ΔCp[i], φ₀[i], Cag[i], t₀[i], τ₀[i], opt; kwargs...)
	end
	return sol
end

"""
    solve_separate(par::Parameters, kwargs...)

Simulate the GC system defined by the structure `par` by solving the
ODEs for ``t(z)`` and ``τ²(z)`` separatly (solving ``t(z)`` and using this result
to solve for ``τ²(z)``) using multiple threads (parallel computing) for the
simulation of different solutes. `kwargs...` to pass additional options to the ODE solve function as named tuples. No multi-threads are used.

Note: The result are two arrays of the solution structure from
DifferentialEquations.jl.

Alternative call:
`solve_separate(L, d, df, gas, T_itp, Fpin_itp, pout_itp, Tchar, θchar, ΔCp, φ₀, Cag, t₀, τ₀, opt; kwargs...)`
with the substance realted quantitites beeing vectors.

# Examples

```julia
julia> sol, peak = solve_separate(par)
```
"""
function solve_separate(par; kwargs...)
    n = length(par.sub)
    sol = Array{Any}(undef, n)
    peak = Array{Any}(undef, n)
    for i=1:n
        sol[i] = solving_migration(par.col, par.prog, par.sub[i], par.opt; kwargs...)
        peak[i] = solving_peakvariance(sol[i], par.col, par.prog, par.sub[i], par.opt; kwargs...)
    end
    return sol, peak
end

function solve_separate(L, d, df, gas, T_itp, Fpin_itp, pout_itp, Tchar, θchar, ΔCp, φ₀, Cag, t₀, τ₀, opt; kwargs...)
    n = length(Tchar)
    sol = Array{Any}(undef, n)
    peak = Array{Any}(undef, n)
    for i=1:n
        sol[i] = solving_migration(L, d, df, gas, T_itp, Fpin_itp, pout_itp, Tchar[i], θchar[i], ΔCp[i], φ₀[i], t₀[i], opt; kwargs...)
        peak[i] = solving_peakvariance(sol[i], L, d, df, gas, T_itp, Fpin_itp, pout_itp, Tchar[i], θchar[i], ΔCp[i], φ₀[i], Cag[i], τ₀[i], opt; kwargs...)
    end
    return sol, peak
end

"""
    solving_migration(col::Column, prog::Program, sub::Substance, opt::Options; kwargs...)

Solve for the migration ``t(z)`` of solute `sub` in the GC Column `col` with
the program `prog` and the options `opt`. `kwargs...` to pass additional options to the ODE solve function as named tuples.

Alternative call:
`solving_migration(L, d, df, gas, T_itp, Fpin_itp, pout_itp, Tchar, θchar, ΔCp, φ₀, t₀, opt; kwargs...)`

Note: The result is the solution structure from
DifferentialEquations.jl.
"""
function solving_migration(col::Column, prog::Program, sub::Substance, opt::Options; kwargs...)
    solution_tz = solving_migration(col.L, col.d, col.df, col.gas, prog.T_itp, prog.Fpin_itp, prog.pout_itp, sub.Tchar, sub.θchar, sub.ΔCp, sub.φ₀, sub.t₀, opt; kwargs...)
    return solution_tz
end

function solving_migration(L, d, df, gas, T_itp, Fpin_itp, pout_itp, Tchar, θchar, ΔCp, φ₀, t₀, opt; kwargs...)
	# this version should be autodiffable for quantities in `p`
    p = (T_itp, Fpin_itp, pout_itp, L, d, df, Tchar, θchar, ΔCp, φ₀, gas, opt) # should be autodiffable for all these parameters
	f_tz(t,p,z) = residency(z, t, p[1], p[2], p[3], p[4], p[5], p[6], p[11], p[7], p[8], p[9], p[10]; ng=p[12].ng, vis=p[12].vis, control=p[12].control, k_th=p[12].k_th)
	t₀ = t₀
	zspan = (0.0, L)
	prob_tz = ODEProblem(f_tz, t₀, zspan, p)
	solution_tz = solve(prob_tz, alg=opt.alg, abstol=opt.abstol, reltol=opt.reltol; kwargs...)
    if SciMLBase.successful_retcode(solution_tz) == false
        @info "ODE solver failed, trying again with set dt=L/1000000."
        solution_tz = solve(prob_tz, alg=opt.alg, abstol=opt.abstol,reltol=opt.reltol; kwargs..., dt=L/1000000)
    end
	return solution_tz
end

"""
    solving_peakvariance(solution_tz, col::Column, prog::Program, sub::Substance, opt::Options; kwargs...)

Solve for the development of the peak variance ``τ²(z)`` of solute `sub` in the GC Column `col` with
the program `prog` and the options `opt` during its migration defined by `solution_tz`. `kwargs...` to pass additional options to the ODE solve function as named tuples.

Alternative call:
`solving_peakvariance(solution_tz, L, d, df, gas, T_itp, Fpin_itp, pout_itp, Tchar, θchar, ΔCp, φ₀, Cag, τ₀, opt; kwargs...)`

Note: The result is the solution structure from
DifferentialEquations.jl.
"""
function solving_peakvariance(solution_tz, col, prog, sub, opt; kwargs...)
    solution_τ²z = solving_peakvariance(solution_tz, col.L, col.d, col.df, col.gas, prog.T_itp, prog.Fpin_itp, prog.pout_itp, sub.Tchar, sub.θchar, sub.ΔCp, sub.φ₀, sub.Cag, sub.τ₀, opt; kwargs...)
    return solution_τ²z
end

function solving_peakvariance(solution_tz, L, d, df, gas, T_itp, Fpin_itp, pout_itp, Tchar, θchar, ΔCp, φ₀, Cag, τ₀, opt; kwargs...)
    # this version should be autodiffable for quantities in `p`
    t(z) = solution_tz(z)
    p = (L, d, df, gas, T_itp, Fpin_itp, pout_itp, Tchar, θchar, ΔCp, φ₀, Cag, opt)
    f_τ²z(τ²,p,z) = peakode(z, t(z), τ², p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[12], p[13])
    τ²₀ = τ₀^2
    zspan = (0.0, L)
    prob_τ²z = ODEProblem(f_τ²z, τ²₀, zspan, p)
    solution_τ²z = solve(prob_τ²z, alg=opt.alg, abstol=opt.abstol,reltol=opt.reltol; kwargs...)
    if SciMLBase.successful_retcode(solution_τ²z) == false
        @info "ODE solver failed, trying again with set dt=L/1000000."
        solution_τ²z = solve(prob_τ²z, alg=opt.alg, abstol=opt.abstol,reltol=opt.reltol; kwargs..., dt=L/1000000)
    end
    return solution_τ²z
end

"""
    solving_odesystem_r(col::Column, prog::Program, sub::Substance, opt::Options; kwargs...)

Solve the migration ``t(z)`` and peak variance development ``τ²(z)`` of solute `sub` in the GC Column `col` with
the program `prog` and the options `opt` as a system of ODEs. `kwargs...` to pass additional options to the ODE solve function as named tuples.

Note: The result is the solution structure from
DifferentialEquations.jl.

See also: [`odesystem_r!`](@ref)
"""
function solving_odesystem_r(col::Column, prog::Program, sub::Substance, opt::Options; kwargs...)
    solution = solving_odesystem_r(col.L, col.d, col.df, col.gas, prog.T_itp, prog.Fpin_itp, prog.pout_itp, sub.Tchar, sub.θchar, sub.ΔCp, sub.φ₀, sub.Cag, sub.t₀, sub.τ₀, opt::GasChromatographySimulator.Options; kwargs...)
    return solution
end

function solving_odesystem_r(L, d, df, gas, T_itp, Fpin_itp, pout_itp, Tchar, θchar, ΔCp, φ₀, Cag, t₀, τ₀, opt::GasChromatographySimulator.Options; kwargs...)
    # this version should be autodiffable for quantities in `p`
    t₀ = [t₀; τ₀^2]
    zspan = (0.0,L)
	p = (L, d, df, gas, T_itp, Fpin_itp, pout_itp, Tchar, θchar, ΔCp, φ₀, Cag, opt)
    prob = ODEProblem(odesystem_r!, t₀, zspan, p)
    solution = solve(prob, alg=opt.alg, abstol=opt.abstol,reltol=opt.reltol; kwargs...)
    if SciMLBase.successful_retcode(solution) == false
        @info "ODE solver failed, trying again with set dt=L/1000000."
        solution = solve(prob, alg=opt.alg, abstol=opt.abstol,reltol=opt.reltol; kwargs..., dt=L/1000000)
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
    L = p[1]
    d = p[2]
    df = p[3]
    gas = p[4]
    T_itp = p[5]
    Fpin_itp = p[6]
    pout_itp = p[7]
    Tchar = p[8]
    θchar = p[9]
    ΔCp = p[10]
    φ₀ = p[11]
    Cag = p[12]
    opt = p[13]
    
    # Compute residency once and reuse it (performance optimization)
    r_val = residency(z, t[1], T_itp, Fpin_itp, pout_itp, L, d, df, gas, Tchar, θchar, ΔCp, φ₀; ng=opt.ng, vis=opt.vis, control=opt.control, k_th=opt.k_th)
    dt[1] = r_val
    # Pass pre-computed r_val to peakode to avoid redundant calculation
    dt[2] = peakode(z, t[1], t[2], L, d, df, gas, T_itp, Fpin_itp, pout_itp, Tchar, θchar, ΔCp, φ₀, Cag, opt, r_val)
end

"""
    peakode(z, t, τ², col, prog, sub, opt)

The second ODE function for the ODE system describing the peak variance
development ``τ²(z)``, using automatic differentiation.

``\\frac{dτ^2}{dz} = H(z, t(z)) r(z, t(z)) + 2 τ^2(z, t(z))
\\frac{∂r}{∂t}(z,t(z))``

See also: [`solving_odesystem_r`](@ref), [`odesystem_r!`](@ref)
"""
function peakode(z, t, τ², col, prog, sub, opt)
    return peakode(z, t, τ², col.L, col.d, col.df, col.gas, prog.T_itp, prog.Fpin_itp, prog.pout_itp, sub.Tchar, sub.θchar, sub.ΔCp, sub.φ₀, sub.Cag, opt)
end

# Optimized version that accepts optional pre-computed residency value to avoid redundant calculations
function peakode(z, t, τ², L, d, df, gas, T_itp, Fpin_itp, pout_itp, Tchar, θchar, ΔCp, φ₀, Cag, opt, r_val=nothing)
    # If r_val is provided, reuse it; otherwise compute it (for backward compatibility)
    if r_val === nothing
        r_val = residency(z, t, T_itp, Fpin_itp, pout_itp, L, d, df, gas, Tchar, θchar, ΔCp, φ₀; ng=opt.ng, vis=opt.vis, control=opt.control, k_th=opt.k_th)
    end
    
    # Compute plate height
    H_val = plate_height(z, t, T_itp, Fpin_itp, pout_itp, L, d, df, gas, Tchar, θchar, ΔCp, φ₀, Cag; ng=opt.ng, vis=opt.vis, control=opt.control, k_th=opt.k_th)
    
    # Compute derivative ∂r/∂t using ForwardDiff
    # ForwardDiff needs to evaluate residency at different points, so we create closures for AD
    if opt.ng==true
        # For ng=true, residency depends on both z and t
        r_func(zt::AbstractVector) = residency(zt[1], zt[2], T_itp, Fpin_itp, pout_itp, L, d, df, gas, Tchar, θchar, ΔCp, φ₀; ng=true, vis=opt.vis, control=opt.control, k_th=opt.k_th)
        ∂r∂t = ForwardDiff.gradient(r_func, [z, t])[2]
    else
        # For ng=false, residency depends only on t (z is fixed)
        r_func_t(t_val) = residency(z, t_val, T_itp, Fpin_itp, pout_itp, L, d, df, gas, Tchar, θchar, ΔCp, φ₀; vis=opt.vis, control=opt.control, k_th=opt.k_th)
        ∂r∂t = ForwardDiff.derivative(r_func_t, t)
    end
    
    return H_val * r_val^2 + 2*τ²*∂r∂t
end
#---End-Solving-Functions---
