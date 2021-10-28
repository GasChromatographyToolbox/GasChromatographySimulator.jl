# make simulation to test this package
using Pkg
Pkg.activate("/Users/janleppert/Documents/GitHub/GasChromatographySimulator")
using GasChromatographySimulator
#using DifferentialEquations
a = 1+ 1
# Options
opt = GasChromatographySimulator.Options(OwrenZen5(), 1e-6, 1e-3, "inlet", true)

# System
a_d = [0.25e-3]
a_df = [0.25e-6]
d(x) = GasChromatographySimulator.gf_const(x, a_d)
df(x) = GasChromatographySimulator.gf_const(x, a_df)
sys = GasChromatographySimulator.System(10.0, d, a_d, df, a_df, "SLB5ms", "He")

# Program
time_steps = [0.0, 60.0, 600.0, 300.0]
temp_steps = [40.0, 40.0, 300.0, 300.0]
pin_steps = [200.0, 200.0, 300.0, 300.0].*1000.0 .+ 101300
pout_steps = [101.3, 101.3, 101.3, 101.3].*1000.0
a_gf = [zeros(length(time_steps)) zeros(length(time_steps)) 10.0.*ones(length(time_steps)) zeros(length(time_steps))]
gf(x) = GasChromatographySimulator.gf_const(x, a_gf)
T_itp = GasChromatographySimulator.temperature_interpolation(time_steps, temp_steps, gf, sys.L)
pin_itp = GasChromatographySimulator.pressure_interpolation(time_steps, pin_steps)
pout_itp = GasChromatographySimulator.pressure_interpolation(time_steps, pout_steps)
prog = GasChromatographySimulator.Program(time_steps, temp_steps, pin_steps, pout_steps, gf, a_gf, T_itp, pin_itp, pout_itp)

# Substance from the load-function-> 1st test
db_path = string(@__DIR__, "/data")
db = "Database_test.csv"
solutes = ["C10", "C11"]
init_t = zeros(length(solutes))
init_τ = zeros(length(solutes))
sub = GasChromatographySimulator.load_solute_database(db_path, db, sys.sp, sys.gas, solutes, init_t, init_τ)
# as test call some values, check the length ...

# Substance as direct input
sub_d = Array{GasChromatographySimulator.Substance}(undef, 2)
sub_d[1] = GasChromatographySimulator.Substance("solute 1",
                                                "0-0-0",
                                                400.0,
                                                30.0,
                                                100.0,
                                                1e-3,
                                                "test",
                                                1e-4,
                                                0.0,
                                                0.0)
sub_d[2] = GasChromatographySimulator.Substance("solute 2",
                                                "1-0-0",
                                                420.0,
                                                31.0,
                                                100.0,
                                                1e-3,
                                                "test",
                                                1e-4,
                                                0.0,
                                                0.0)

# parameters
par = GasChromatographySimulator.Par(sys, prog, sub, opt)
par_d = GasChromatographySimulator.Par(sys, prog, sub_d, opt)

# solving
# for odesys:
    sol_odesys = GasChromatographySimulator.solve_system_multithreads(par)
    pl_odesys = GasChromatographySimulator.peaklist(sol_odesys, par)
    # test: tR(C10) = 97.943, τR(C10) = 0.581
    # test: tR(C11) = 140.416, τR(C11) = 0.540
    round.(pl_odesys.tR, digits=3) == [97.943, 140.416]
    round.(pl_odesys.τR, digits=3) == [0.581, 0.540]
    
    sol_odesys_ng = GasChromatographySimulator.solve_system_multithreads(par, ng=true)
    pl_odesys_ng = GasChromatographySimulator.peaklist(sol_odesys_ng, par)
    # non-gradient separation should give the 'same' results with both methods (ng=false and ng=true)
    abs.(1 .- pl_odesys_ng.tR./pl_odesys.tR)[1] < 1e-5
    abs.(1 .- pl_odesys_ng.tR./pl_odesys.tR)[2] < 1e-5
    abs.(1 .- pl_odesys_ng.τR./pl_odesys.τR)[1] < 1e-4
    abs.(1 .- pl_odesys_ng.τR./pl_odesys.τR)[2] < 1e-4

    sol_odesys_d = GasChromatographySimulator.solve_system_multithreads(par_d)
    pl_odesys_d = GasChromatographySimulator.peaklist(sol_odesys_d, par_d)
    # test: tR(C10) = 145.656, τR^2(C10) = 0.506
    # test: tR(C11) = 188.478, τR^2(C11) = 0.292
    round.(pl_odesys_d.tR, digits=3) == [145.656, 188.478]
    round.(pl_odesys_d.τR, digits=3) == [0.506, 0.483]
    
    sol_odesys_d_ng = GasChromatographySimulator.solve_system_multithreads(par_d, ng=true)
    pl_odesys_d_ng = GasChromatographySimulator.peaklist(sol_odesys_d_ng, par_d)
    # non-gradient separation should give the 'same' results with both methods (ng=false and ng=true)
    abs.(1 .- pl_odesys_d_ng.tR./pl_odesys_d.tR)[1] < 1e-5
    abs.(1 .- pl_odesys_d_ng.tR./pl_odesys_d.tR)[2] < 1e-5
    abs.(1 .- pl_odesys_d_ng.τR./pl_odesys_d.τR)[1] < 1e-4
    abs.(1 .- pl_odesys_d_ng.τR./pl_odesys_d.τR)[2] < 1e-4

# for migration and peakvariance separatly:
    sol_migr, sol_peak = GasChromatographySimulator.solve_multithreads(par)
    pl_migr = GasChromatographySimulator.peaklist(sol_migr, sol_peak, par)
    # test: tR(C10) = 97.944, τR(C10) = 0.578
    # test: tR(C11) = 140.403, τR(C11) = 0.540
    round.(pl_migr.tR, digits=3) == [97.944, 140.403]
    round.(pl_migr.τR, digits=3) == [0.578, 0.540]

    sol_migr_ng, sol_peak_ng = GasChromatographySimulator.solve_multithreads(par, ng=true)
    pl_migr_ng = GasChromatographySimulator.peaklist(sol_migr_ng, sol_peak_ng, par)
    abs.(1 .- pl_migr_ng.tR./pl_migr.tR)[1] < 1e-5
    abs.(1 .- pl_migr_ng.tR./pl_migr.tR)[2] < 1e-5
    abs.(1 .- pl_migr_ng.τR./pl_migr.τR)[1] < 1e-4
    abs.(1 .- pl_migr_ng.τR./pl_migr.τR)[2] < 1e-4