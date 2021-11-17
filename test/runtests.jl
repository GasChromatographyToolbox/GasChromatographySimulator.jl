using Test, GasChromatographySimulator

# Options
opt = GasChromatographySimulator.Options(OwrenZen5(), 1e-6, 1e-3, "inlet", true)

# System
L = 10.0
a_d = [0.25e-3]
a_df = [0.25e-6]
d(x) = GasChromatographySimulator.gradient(x, a_d)
df(x) = GasChromatographySimulator.gradient(x, a_df)
sp = "SLB5ms"
gas = "He"
sys = GasChromatographySimulator.System(10.0, d, a_d, df, a_df, sp, gas)

sys_c = GasChromatographySimulator.constructor_System(L, a_d[1], a_df[1], sp, gas)

@test sys.d(sys.L) == sys_c.d(sys_c.L)

# Program
time_steps = [0.0, 60.0, 600.0, 300.0]
temp_steps = [40.0, 40.0, 300.0, 300.0]
pin_steps = [200.0, 200.0, 300.0, 300.0].*1000.0 .+ 101300
pout_steps = [101.3, 101.3, 101.3, 101.3].*1000.0
a_gf = [zeros(length(time_steps)) zeros(length(time_steps)) sys.L.*ones(length(time_steps)) zeros(length(time_steps))]
gf(x) = GasChromatographySimulator.gradient(x, a_gf)
T_itp = GasChromatographySimulator.temperature_interpolation(time_steps, temp_steps, gf, sys.L)
pin_itp = GasChromatographySimulator.pressure_interpolation(time_steps, pin_steps)
pout_itp = GasChromatographySimulator.pressure_interpolation(time_steps, pout_steps)
prog = GasChromatographySimulator.Program(time_steps, temp_steps, pin_steps, pout_steps, gf, a_gf, T_itp, pin_itp, pout_itp)

prog_c = GasChromatographySimulator.constructor_Program(time_steps, temp_steps, pin_steps, pout_steps, ΔT_steps, x₀_steps, L₀_steps, α_steps, opt.Tcontrol, sys.L)
prog_c_ng = GasChromatographySimulator.constructor_Program(time_steps, temp_steps, pin_steps, pout_steps, sys.L)

@test prog.T_itp == prog_c.T_itp
@test prog_c.T_itp == prog_c_ng.T_itp

# Substance from the load-function-> 1st test
db_path = string(@__DIR__, "/data")
db = "Database_test.csv"
solutes = ["C10", "C11"]
init_t = zeros(length(solutes))
init_τ = zeros(length(solutes))
sub = GasChromatographySimulator.load_solute_database(db_path, db, sys.sp, sys.gas, solutes, init_t, init_τ)

# parameters
par = GasChromatographySimulator.Parameters(sys, prog, sub, opt)

# solving
# for odesys:
sol_odesys = GasChromatographySimulator.solve_system_multithreads(par)
pl_odesys = GasChromatographySimulator.peaklist(sol_odesys, par)
# test: tR(C10) = 97.943, τR(C10) = 0.581
# test: tR(C11) = 140.416, τR(C11) = 0.540
@test round.(pl_odesys.tR, digits=3) == [97.943, 140.416]
@test round.(pl_odesys.τR, digits=3) == [0.581, 0.540]

sol_odesys_ng = GasChromatographySimulator.solve_system_multithreads(par, ng=true)
pl_odesys_ng = GasChromatographySimulator.peaklist(sol_odesys_ng, par)
# non-gradient separation should give the 'same' results with both methods (ng=false and ng=true)
@test abs.(1 .- pl_odesys_ng.tR./pl_odesys.tR)[1] < 1e-5
@test abs.(1 .- pl_odesys_ng.τR./pl_odesys.τR)[2] < 1e-4

# for migration and peakvariance separatly:
sol_migr, sol_peak = GasChromatographySimulator.solve_multithreads(par)
pl_migr = GasChromatographySimulator.peaklist(sol_migr, sol_peak, par)
# test: tR(C10) = 97.944, τR(C10) = 0.578
# test: tR(C11) = 140.403, τR(C11) = 0.540
@test round.(pl_migr.tR, digits=3) == [97.944, 140.403]
@test round.(pl_migr.τR, digits=3) == [0.578, 0.540]

sol_migr_ng, sol_peak_ng = GasChromatographySimulator.solve_multithreads(par, ng=true)
pl_migr_ng = GasChromatographySimulator.peaklist(sol_migr_ng, sol_peak_ng, par)
@test abs.(1 .- pl_migr_ng.tR./pl_migr.tR)[1] < 1e-5
@test abs.(1 .- pl_migr_ng.τR./pl_migr.τR)[2] < 1e-4

# test of holdup_time:
t = rand()*cumsum(time_steps)[end]
T_test = par.prog.T_itp(0.0, t)
η_T = GasChromatographySimulator.viscosity(T_test, "He")
η_t = GasChromatographySimulator.viscosity(0.0, t, par.prog.T_itp, par.sys.gas)

@test η_T == η_t

tM_T = GasChromatographySimulator.holdup_time(T_test, par.prog.pin_itp(t), par.prog.pout_itp(t), L, a_d[1], gas)
tM_t = GasChromatographySimulator.holdup_time(t, par.prog.T_itp, par.prog.pin_itp, par.prog.pout_itp, par.sys.L, par.sys.d, par.sys.gas; ng=false)
tM_t_ng = GasChromatographySimulator.holdup_time(t, par.prog.T_itp, par.prog.pin_itp, par.prog.pout_itp, par.sys.L, par.sys.d, par.sys.gas; ng=true)

@test (tM_T - tM_t)/tM_t < 1e-10
@test (tM_t_ng - tM_t)/tM_t < 1e-10