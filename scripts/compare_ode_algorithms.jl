#!/usr/bin/env julia
#=
Comparison script for ODE solver algorithms
Tests multiple algorithms using Database_test.csv substances
Includes two test cases:
1. Without thermal gradient (long column)
2. With thermal gradient (short column)

Algorithms tested:
- OwrenZen5() (reference/baseline)
- Tsit5() (modern default, recommended)
- Vern9() (high accuracy)
- BS5() (Bogacki-Shampine 5)
- OwrenZen4() (faster alternative)
- DP5() (Dormand-Prince 5, classic)
=#

using GasChromatographySimulator
using OrdinaryDiffEq: Tsit5, Vern9, BS5, OwrenZen4, DP5
using Printf
using Statistics  # For mean() function

# Use BenchmarkTools for detailed benchmarking
using BenchmarkTools
HAS_BENCHMARKTOOLS = true

# Define algorithms to test (reference algorithm should be first)
ALGORITHMS = [
    ("OwrenZen5", OwrenZen5()),  # Reference/baseline
    ("Tsit5", Tsit5()),
    ("Vern9", Vern9()),
    ("BS5", BS5()),
    ("OwrenZen4", OwrenZen4()),
    ("DP5", DP5()),
]

# Configuration
db_path = abspath(joinpath(@__DIR__, "..", "data"))
db_name = "Database_test.csv"
stat_phase = "Wax"  # Stationary phase to use

# Load database to get available solutes
using CSV, DataFrames
db_dataframe = DataFrame(CSV.File(joinpath(db_path, db_name), header=1, silencewarnings=true, stringtype=String))
solutes = GasChromatographySimulator.all_solutes(stat_phase, db_dataframe)

println("="^80)
println("ODE Solver Algorithm Comparison")
println("Testing $(length(ALGORITHMS)) algorithms: ", join([a[1] for a in ALGORITHMS], ", "))
println("="^80)
println()
println("Loaded $(length(solutes)) substances from database for phase: $stat_phase")
println("Substances: ", join(solutes[1:min(5, length(solutes))], ", "), length(solutes) > 5 ? "..." : "")
println()

# Load substances
t₀ = zeros(length(solutes))
τ₀ = zeros(length(solutes))
sub = GasChromatographySimulator.load_solute_database(db_path, db_name, stat_phase, "He", solutes, t₀, τ₀)

# Function to run comparison for a test case with multiple algorithms
function run_comparison(case_name, col, prog, algorithms, sub, solutes)
    println("="^80)
    println(case_name)
    println("="^80)
    println()
    
    println("Column: L=$(col.L)m, d=$(col.d*1000)mm, df=$(col.df*1e6)μm, gas=$(col.gas)")
    has_gradient = hasfield(typeof(prog), :a_gf) && !isempty(prog.a_gf) && size(prog.a_gf, 1) > 0 && any(abs.(prog.a_gf[:, 1]) .> 1e-10)
    if has_gradient
        println("Program: $(length(prog.time_steps)) steps, $(prog.temp_steps[1])°C → $(prog.temp_steps[end])°C")
        println("Thermal gradient: ΔT=$(prog.a_gf[1,1])°C (initial), α=$(prog.a_gf[1,4]) → $(prog.a_gf[end,4])")
    else
        println("Program: $(length(prog.time_steps)) steps, $(prog.temp_steps[1])°C → $(prog.temp_steps[end])°C (no gradient)")
    end
    println()
    
    # Create options and parameters for all algorithms
    opts = Dict()
    pars = Dict()
    for (name, alg) in algorithms
        opts[name] = GasChromatographySimulator.Options(alg=alg, abstol=1e-6, reltol=1e-3)
        pars[name] = GasChromatographySimulator.Parameters(col, prog, sub, opts[name])
    end
    
    # Reference algorithm (first one)
    ref_name = algorithms[1][1]
    
    # Run simulations for all algorithms
    println("Running simulations...")
    println()
    
    solutions = Dict()
    peaklists = Dict()
    times = Dict()
    
    for (name, _) in algorithms
        println("Running $(name)() simulation...")
        t_start = time()
        solutions[name] = GasChromatographySimulator.solve_system_multithreads(pars[name])
        times[name] = time() - t_start
        peaklists[name] = GasChromatographySimulator.peaklist(solutions[name], pars[name])
        println("  Completed in $(@sprintf("%.3f", times[name])) s")
    end
    println()
    
    # Check for successful solutions
    println("Solution Status:")
    success_flags = Dict()
    for (name, _) in algorithms
        success_flags[name] = all(s -> GasChromatographySimulator.SciMLBase.successful_retcode(s), solutions[name])
        status = success_flags[name] ? "✓ All successful" : "✗ Some failed"
        println(@sprintf("  %-15s: %s", name*"()", status))
    end
    println()
    
    # Compare retention times against reference
    println("Retention Time Comparison (tR in seconds, differences vs $(ref_name)()):")
    println("-"^100)
    header = @sprintf("%-25s", "Substance")
    for (name, _) in algorithms
        header *= @sprintf(" %15s", name)
    end
    header *= @sprintf(" %15s", "Max Diff")
    println(header)
    println("-"^100)
    
    max_diffs_tR = Dict()
    max_diff_idx_tR = Dict()
    for (name, _) in algorithms
        if name != ref_name
            max_diffs_tR[name] = 0.0
            max_diff_idx_tR[name] = 0
        end
    end
    
    for i in 1:length(solutes)
        ref_tR = peaklists[ref_name].tR[i]
        if !ismissing(ref_tR)
            row = @sprintf("%-25s", solutes[i])
            for (name, _) in algorithms
                tR_val = peaklists[name].tR[i]
                if !ismissing(tR_val)
                    row *= @sprintf(" %15.6f", tR_val)
                    if name != ref_name
                        diff = abs(tR_val - ref_tR)
                        if diff > max_diffs_tR[name]
                            max_diffs_tR[name] = diff
                            max_diff_idx_tR[name] = i
                        end
                    end
                else
                    row *= @sprintf(" %15s", "Failed")
                end
            end
            # Calculate max difference across all algorithms
            max_diff_all = 0.0
            for (name, _) in algorithms
                if name != ref_name && !ismissing(peaklists[name].tR[i])
                    diff = abs(peaklists[name].tR[i] - ref_tR)
                    max_diff_all = max(max_diff_all, diff)
                end
            end
            row *= @sprintf(" %15.6e", max_diff_all)
            println(row)
        end
    end
    
    println("-"^100)
    println("Maximum differences vs $(ref_name)():")
    for (name, _) in algorithms
        if name != ref_name && haskey(max_diffs_tR, name) && max_diff_idx_tR[name] > 0
            println(@sprintf("  %-15s: %.6e s (substance: %s)", name*"()", max_diffs_tR[name], solutes[max_diff_idx_tR[name]]))
        end
    end
    println()
    
    # Performance comparison
    println("="^80)
    println("PERFORMANCE COMPARISON")
    println("="^80)
    println()

    bench_median = Dict{String, Float64}()
    bench_min = Dict{String, Float64}()
    bench_max = Dict{String, Float64}()

    if HAS_BENCHMARKTOOLS
        println("Benchmarking algorithms (3 samples each)...")
        benchmarks = Dict()
        for (name, _) in algorithms
            println("  $(name)()...")
            benchmarks[name] = @benchmark GasChromatographySimulator.solve_system_multithreads($(pars[name])) samples=3 evals=1
        end
        
        println()
        println("Performance Results (median time):")
        println("-"^80)
        println(@sprintf("%-15s %15s %15s %15s %15s", "Algorithm", "Time (s)", "vs Ref", "Speedup", "Range (s)"))
        println("-"^80)
        
        ref_time = median(benchmarks[ref_name].times) / 1e9
        for (name, _) in algorithms
            med_time = median(benchmarks[name].times) / 1e9
            min_time = minimum(benchmarks[name].times) / 1e9
            max_time = maximum(benchmarks[name].times) / 1e9
            bench_median[name] = med_time
            bench_min[name] = min_time
            bench_max[name] = max_time
            speedup = ref_time / med_time
            speedup_str = name == ref_name ? "1.00x" : @sprintf("%.2fx", speedup)
            vs_ref = name == ref_name ? "ref" : @sprintf("%+.1f%%", (med_time - ref_time) / ref_time * 100)
            println(@sprintf("%-15s %15.3f %15s %15s %15s", name*"()", med_time, vs_ref, speedup_str, "$(@sprintf("%.3f", min_time))-$(@sprintf("%.3f", max_time))"))
        end
        println("-"^80)
    else
        println("Performance Results (simple timing):")
        println("-"^80)
        println(@sprintf("%-15s %15s", "Algorithm", "Time (s)"))
        println("-"^80)
        for (name, _) in algorithms
            println(@sprintf("%-15s %15.3f", name*"()", times[name]))
            bench_median[name] = times[name]
            bench_min[name] = times[name]
            bench_max[name] = times[name]
        end
        println("-"^80)
    end
    println()
    
    # ODE solver statistics
    println("="^80)
    println("ODE SOLVER STATISTICS")
    println("="^80)
    println()
    
    println(@sprintf("%-15s %15s %15s %15s %15s", "Algorithm", "Avg Steps", "Avg FEvals", "Total Steps", "Total FEvals"))
    println("-"^80)
    
    for (name, _) in algorithms
        steps = [length(s.t) for s in solutions[name]]
        fevals = [s.stats.nf for s in solutions[name]]
        println(@sprintf("%-15s %15.1f %15.1f %15d %15d", 
            name*"()", mean(steps), mean(fevals), sum(steps), sum(fevals)))
    end
    println()
    
    # Summary
    println("="^80)
    println("SUMMARY")
    println("="^80)
    println()
    
    println("Numerical Accuracy (max difference vs $(ref_name)()):")
    for (name, _) in algorithms
        if name != ref_name && haskey(max_diffs_tR, name)
            diff = max_diffs_tR[name]
            if diff < 1e-6
                println("  $(name)(): ✓ Excellent (< 1e-6 s)")
            elseif diff < 1e-3
                println("  $(name)(): ✓ Good (< 1e-3 s)")
            else
                println("  $(name)(): ⚠ Significant (> 1e-3 s): $(@sprintf("%.6e", diff)) s")
            end
        end
    end
    println()
    
    println("Reliability:")
    for (name, _) in algorithms
        status = success_flags[name] ? "✓ All successful" : "✗ Some failures"
        println("  $(name)(): $status")
    end
    println()
    
    # Return results
    results = Dict()
    for (name, _) in algorithms
        results[name] = (
            max_diff_tR = name == ref_name ? 0.0 : (haskey(max_diffs_tR, name) ? max_diffs_tR[name] : Inf),
            success = success_flags[name],
            solution = solutions[name],
            peaklist = peaklists[name],
            time = haskey(times, name) ? times[name] : 0.0,
            bench_median_s = haskey(bench_median, name) ? bench_median[name] : NaN,
            bench_min_s = haskey(bench_min, name) ? bench_min[name] : NaN,
            bench_max_s = haskey(bench_max, name) ? bench_max[name] : NaN
        )
    end
    
    return results
end

# ============================================================================
# TEST CASE 1: Without thermal gradient (long column)
# ============================================================================
L1 = 30.0  # Column length in m
d1 = 0.25e-3  # Column diameter in m
df1 = 0.25e-6  # Film thickness in m
gas1 = "He"  # Carrier gas

time_steps_1 = [0.0, 60.0, 1800.0, 120.0]  # Hold at 40°C for 1 min, then ramp to 300°C
temp_steps_1 = [40.0, 40.0, 300.0, 300.0]  # Temperature in °C
Fpin_steps_1 = [300000.0, 300000.0, 450000.0, 450000.0]  # Inlet pressure in Pa
pout_steps_1 = [101300.0, 101300.0, 101300.0, 101300.0]  # Outlet pressure in Pa

col1 = GasChromatographySimulator.Column(L1, d1, df1, stat_phase, gas1)
prog1 = GasChromatographySimulator.Program(time_steps_1, temp_steps_1, Fpin_steps_1, pout_steps_1, L1)

results1 = run_comparison("TEST CASE 1: Without Thermal Gradient (Long Column, L=$(L1)m)", 
                          col1, prog1, ALGORITHMS, sub, solutes)

# ============================================================================
# TEST CASE 2: With thermal gradient (short column)
# ============================================================================
L2 = 2.0  # Short column length in m
d2 = 0.1e-3  # Column diameter in m
df2 = 0.1e-6  # Film thickness in m
gas2 = "He"  # Carrier gas

# Temperature program with thermal gradient
time_steps_2 = [0.0, 60.0, 300.0, 120.0]  # Time steps in s
temp_steps_2 = [40.0, 40.0, 320.0, 320.0]  # Temperature steps in °C
Fpin_steps_2 = [300000.0, 300000.0, 450000.0, 450000.0]  # Inlet pressure in Pa
pout_steps_2 = [101300.0, 101300.0, 101300.0, 101300.0]  # Outlet pressure in Pa
ΔT_steps_2 = [50.0, 50.0, 50.0, 50.0]  # Temperature difference in °C
x₀_steps_2 = [0.0, 0.0, 0.0, 0.0]  # Spatial offset of gradient in m
L₀_steps_2 = [L2, L2, L2, L2]  # Distance over which ΔT is measured in m
α_steps_2 = [0.0, -2.0, -5.0, -5.0]  # Gradient profile factor
Tcontrol_2 = "inlet"  # Temperature control point

col2 = GasChromatographySimulator.Column(L2, d2, df2, stat_phase, gas2)
prog2 = GasChromatographySimulator.Program(time_steps_2, temp_steps_2, Fpin_steps_2, pout_steps_2, 
                                           ΔT_steps_2, x₀_steps_2, L₀_steps_2, α_steps_2, Tcontrol_2, L2)

results2 = run_comparison("TEST CASE 2: With Thermal Gradient (Short Column, L=$(L2)m)", 
                          col2, prog2, ALGORITHMS, sub, solutes)

# ============================================================================
# OVERALL SUMMARY AND RECOMMENDATIONS
# ============================================================================
println("="^80)
println("OVERALL SUMMARY AND RECOMMENDATIONS")
println("="^80)
println()

ref_name = ALGORITHMS[1][1]

println("Test Case 1 (No Gradient) - Max differences vs $(ref_name)():")
for (name, _) in ALGORITHMS
    if name != ref_name && haskey(results1, name)
        println(@sprintf("  %-15s: %.6e s, %s", name*"()", results1[name].max_diff_tR, 
            results1[name].success ? "✓" : "✗"))
    end
end
println()

println("Test Case 2 (With Gradient) - Max differences vs $(ref_name)():")
for (name, _) in ALGORITHMS
    if name != ref_name && haskey(results2, name)
        println(@sprintf("  %-15s: %.6e s, %s", name*"()", results2[name].max_diff_tR, 
            results2[name].success ? "✓" : "✗"))
    end
end
println()

println("Overall Assessment:")
# Find best performers
best_accuracy = Inf
best_accuracy_name = ""
best_speed = 0.0
best_speed_name = ""

for (name, _) in ALGORITHMS
    if name != ref_name && haskey(results1, name) && haskey(results2, name)
        # Combined accuracy (max of both test cases)
        combined_diff = max(results1[name].max_diff_tR, results2[name].max_diff_tR)
        if combined_diff < best_accuracy && results1[name].success && results2[name].success
            global best_accuracy = combined_diff
            global best_accuracy_name = name
        end
        
        # Speed (using test case 1)
        if haskey(results1[name], :time) && results1[name].time > 0
            speed = 1.0 / results1[name].time
            if speed > best_speed && results1[name].success
                global best_speed = speed
                global best_speed_name = name
            end
        end
    end
end

if best_accuracy < Inf
    println("  Best accuracy: $(best_accuracy_name)() (max diff: $(@sprintf("%.6e", best_accuracy)) s)")
end
if best_speed_name != ""
    println("  Fastest: $(best_speed_name)()")
end

# Check if any algorithm shows significantly worse accuracy in gradient case
println()
println("Gradient Case Analysis:")
for (name, _) in ALGORITHMS
    if name != ref_name && haskey(results1, name) && haskey(results2, name)
        ratio = results2[name].max_diff_tR / (results1[name].max_diff_tR + 1e-10)
        if ratio > 10.0
            println("  ⚠ $(name)() shows $(@sprintf("%.1f", ratio))x larger differences in gradient case")
        elseif ratio < 0.1
            println("  ✓ $(name)() performs better in gradient case")
        end
    end
end
println()

# Optional export of results for documentation plots
out_path = get(ENV, "GCS_BENCH_OUT", "")
if !isempty(out_path)
    function _build_results_df(case_label, results, algorithms, ref_name)
        rows = DataFrame(
            case = String[],
            algorithm = String[],
            ref_algorithm = String[],
            max_diff_tR_s = Float64[],
            success = Bool[],
            run_time_s = Float64[],
            bench_median_s = Float64[],
            bench_min_s = Float64[],
            bench_max_s = Float64[],
        )
        for (name, _) in algorithms
            r = results[name]
            push!(rows, (
                case_label,
                name,
                ref_name,
                r.max_diff_tR,
                r.success,
                r.time,
                r.bench_median_s,
                r.bench_min_s,
                r.bench_max_s,
            ))
        end
        return rows
    end

    mkpath(dirname(out_path))
    df_out = vcat(
        _build_results_df("no_gradient", results1, ALGORITHMS, ref_name),
        _build_results_df("gradient", results2, ALGORITHMS, ref_name),
    )
    CSV.write(out_path, df_out)
    println("Wrote benchmark results to: $(out_path)")
end

println("="^80)
println("Comparison complete!")
println("="^80)
