#!/usr/bin/env julia
#=
Speed comparison of the time-domain solvers and the spatial-domain solvers
on the same forward elution.

Both ends have a finite pressure. A vacuum outlet cannot reverse, so it is
not an equivalent case for the time-domain solver.

Pairs, each timed as one call:
  migration     solving_migration              vs solving_migration_time
  separate      migration + peakvariance       vs migration_time + bandvariance_time
  combined      solving_odesystem_r            vs solving_odesystem_time

The time-domain functions are called as a user would call them: default
t_end (last program knot), exit callback, and program knots as tstops.
The spatial functions have no callback and no tstops.

  julia --project scripts/benchmark_time_vs_spatial.jl
=#

using GasChromatographySimulator
using BenchmarkTools
using Printf
using Statistics

const GCS = GasChromatographySimulator

struct Case
    name::String
    col::GCS.Column
    prog::GCS.Program
    sub::GCS.Substance
    opt::GCS.Options
end

function isothermal_case()
    L = 1.0
    col = GCS.Column(L, 0.25e-3, 0.25e-6, "DB5", "He")
    prog = GCS.Program([0.0, 60.0, 120.0], fill(100.0, 3), fill(150000.0, 3), fill(120000.0, 3))
    sub = GCS.Substance("test", "1-1-1", 400.0, 30.0, 100.0, 1e-3, "bench", 1e-6, 0.0, 0.0)
    return Case("isothermal 1 m, 100 °C, 150/120 kPa", col, prog, sub, GCS.Options())
end

function programmed_case()
    L = 10.0
    col = GCS.Column(L, 0.25e-3, 0.25e-6, "DB5", "He")
    # Durations. Absolute knots are 0, 60, 660, 780 s. Outlet is atmospheric.
    prog = GCS.Program(
        [0.0, 60.0, 600.0, 120.0],
        [40.0, 40.0, 300.0, 300.0],
        [18.0, 18.0, 98.0, 98.0] .* 1000.0 .+ 101300.0,
        fill(101300.0, 4),
    )
    sub = GCS.Substance("test", "1-1-1", 400.0, 30.0, 100.0, 1e-3, "bench", 1e-6, 0.0, 0.0)
    return Case("programmed 10 m, 40→300 °C, outlet 101.3 kPa", col, prog, sub, GCS.Options())
end

nf(sol) = sol.destats.nf
naccept(sol) = sol.destats.naccept
nreject(sol) = sol.destats.nreject

function retention(sol, domain)
    if domain == :spatial
        return sol.u[end] isa Number ? sol.u[end] : sol.u[end][1]
    else
        return sol.t[end]
    end
end

function temporal_width(sol_pos, sol_time, case)
    τ_pos = sqrt(sol_pos.u[end] isa Number ? sol_pos.u[end] : sol_pos.u[end][2])
    σ² = sol_time.u[end] isa Number ? sol_time.u[end] : sol_time.u[end][2]
    u_exit = GCS.solute_velocity(case.col.L, sol_time.t[end], case.col, case.prog, case.sub, case.opt)
    return τ_pos, sqrt(σ²) / abs(u_exit)
end

function summarize(label, b, evals, accepted, rejected)
    med = median(b.times) / 1e6
    mn = minimum(b.times) / 1e6
    allocs = b.allocs
    mem = b.memory / 1024
    per = evals == 0 ? NaN : med * 1e3 / evals
    @printf("  %-22s %8.3f ms   min %8.3f   %6d evals   %5d accept   %4d reject   %7.2f µs/eval   %6d allocs   %7.1f KiB\n",
        label, med, mn, evals, accepted, rejected, per, allocs, mem)
    return med
end

function run_case(case)
    println()
    println("=" ^ 100)
    println(case.name)
    println("Options: $(typeof(case.opt.alg).name.name), abstol=$(case.opt.abstol), reltol=$(case.opt.reltol), ng=$(case.opt.ng)")
    knots = collect(GCS.program_knots(case.prog.T_itp, case.sub.t₀, case.prog.T_itp.grid[2][end]))
    println("Program knots inside the default time span: ", isempty(knots) ? "none" : join(knots, ", "))
    println("-" ^ 100)

    col, prog, sub, opt = case.col, case.prog, case.sub, case.opt

    spatial_mig() = GCS.solving_migration(col, prog, sub, opt)
    time_mig() = GCS.solving_migration_time(col, prog, sub, opt)
    function spatial_sep()
        sol = GCS.solving_migration(col, prog, sub, opt)
        var = GCS.solving_peakvariance(sol, col, prog, sub, opt)
        return sol, var
    end
    function time_sep()
        sol = GCS.solving_migration_time(col, prog, sub, opt)
        var = GCS.solving_bandvariance_time(sol, col, prog, sub, opt)
        return sol, var
    end
    spatial_sys() = GCS.solving_odesystem_r(col, prog, sub, opt)
    time_sys() = GCS.solving_odesystem_time(col, prog, sub, opt)

    # Compile, then check that both domains elute at the same time.
    sm, tm = spatial_mig(), time_mig()
    ss, ts = spatial_sep(), time_sep()
    sr, tr = spatial_sys(), time_sys()

    tR_s, tR_t = retention(sm, :spatial), retention(tm, :time)
    println(@sprintf("Migration check:    spatial t(L) = %.6f s   time-domain t_exit = %.6f s   |Δ| = %.3e s   x_end = %.6f m   retcodes %s / %s",
        tR_s, tR_t, abs(tR_s - tR_t), tm.u[end], sm.retcode, tm.retcode))
    knots_before_exit = knots[knots .< tR_t]
    println("Knots the time-domain solver must hit before exit: ", isempty(knots_before_exit) ? "none" : join(knots_before_exit, ", "))
    τs, τt = temporal_width(ss[2], ts[2], case)
    println(@sprintf("Separate variance:  spatial τ = %.6e s   time-domain τ = %.6e s   |Δ| = %.3e s", τs, τt, abs(τs - τt)))
    τsr, τtr = temporal_width(sr, tr, case)
    println(@sprintf("Combined system:    spatial τ = %.6e s   time-domain τ = %.6e s   |Δ| = %.3e s   retcodes %s / %s",
        τsr, τtr, abs(τsr - τtr), sr.retcode, tr.retcode))
    println()

    samples = 10
    seconds = 15
    println(@sprintf("%-24s %s", "", "median"))

    b = @benchmark $(spatial_mig)() samples=samples seconds=seconds evals=1
    med_s = summarize("spatial migration", b, nf(sm), naccept(sm), nreject(sm))
    b = @benchmark $(time_mig)() samples=samples seconds=seconds evals=1
    med_t = summarize("time migration", b, nf(tm), naccept(tm), nreject(tm))
    @printf("  ratio time/spatial (migration)            %.2f\n", med_t / med_s)

    b = @benchmark $(spatial_sep)() samples=samples seconds=seconds evals=1
    med_s = summarize("spatial separate", b, nf(ss[1]) + nf(ss[2]), naccept(ss[1]) + naccept(ss[2]), nreject(ss[1]) + nreject(ss[2]))
    b = @benchmark $(time_sep)() samples=samples seconds=seconds evals=1
    med_t = summarize("time separate", b, nf(ts[1]) + nf(ts[2]), naccept(ts[1]) + naccept(ts[2]), nreject(ts[1]) + nreject(ts[2]))
    @printf("  ratio time/spatial (separate)             %.2f\n", med_t / med_s)

    b = @benchmark $(spatial_sys)() samples=samples seconds=seconds evals=1
    med_s = summarize("spatial combined", b, nf(sr), naccept(sr), nreject(sr))
    b = @benchmark $(time_sys)() samples=samples seconds=seconds evals=1
    med_t = summarize("time combined", b, nf(tr), naccept(tr), nreject(tr))
    @printf("  ratio time/spatial (combined)             %.2f\n", med_t / med_s)
end

println("Time-domain vs spatial-domain solver benchmark")
println("Julia $(VERSION), GasChromatographySimulator $(pkgversion(GCS))")
run_case(isothermal_case())
run_case(programmed_case())
println()
println("Done.")
