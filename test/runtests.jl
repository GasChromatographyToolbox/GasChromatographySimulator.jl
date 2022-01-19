using Test, GasChromatographySimulator

# general parameters
L = 10.0
d = 0.25e-3
df = 0.25e-6
sp = "SLB5ms"
gas = "He"
time_steps = [0.0, 60.0, 600.0, 300.0]
temp_steps = [40.0, 40.0, 300.0, 300.0]
pin_steps = [200.0, 200.0, 300.0, 300.0].*1000.0 .+ 101300
pout_steps = [0.0, 0.0, 0.0, 0.0].*1000.0
ΔT_steps = [20.0, 30.0, 50.0, 40.0]
x₀_steps = zeros(length(time_steps))
L₀_steps = L.*ones(length(time_steps))
α_steps = zeros(length(time_steps))
db_path = string(@__DIR__, "/data")
db = "Database_test.csv"
@testset "parameters check" begin
    # default Options
    opt = GasChromatographySimulator.Options(OwrenZen5(), 1e-6, 1e-3, "inlet", true)
    @test GasChromatographySimulator.Options() == opt

    # System
    a_d = [d]
    a_df = [df]
    d_fun(x) = GasChromatographySimulator.gradient(x, a_d)
    df_fun(x) = GasChromatographySimulator.gradient(x, a_df)
    sys = GasChromatographySimulator.System(L, d_fun, a_d, df_fun, a_df, sp, gas)
    sys_c = GasChromatographySimulator.constructor_System(L, a_d[1], a_df[1], sp, gas)
    @test sys.d(sys.L) == sys_c.d(sys_c.L)

    # Program 
    a_gf = [ΔT_steps x₀_steps L₀_steps α_steps]
    gf(x) = GasChromatographySimulator.gradient(x, a_gf)
    T_itp = GasChromatographySimulator.temperature_interpolation(time_steps, temp_steps, gf, sys.L)
    pin_itp = GasChromatographySimulator.pressure_interpolation(time_steps, pin_steps)
    pout_itp = GasChromatographySimulator.pressure_interpolation(time_steps, pout_steps)
    prog = GasChromatographySimulator.Program(time_steps, temp_steps, pin_steps, pout_steps, gf, a_gf, T_itp, pin_itp, pout_itp)
    prog_c = GasChromatographySimulator.constructor_Program(time_steps, temp_steps, pin_steps, pout_steps, ΔT_steps, x₀_steps, L₀_steps, α_steps, opt.Tcontrol, sys.L)
    prog_c_ng = GasChromatographySimulator.constructor_Program(time_steps, temp_steps, pin_steps, pout_steps, sys.L) # without gradient
    @test prog.T_itp == prog_c.T_itp
    @test prog_c.T_itp != prog_c_ng.T_itp
    @test prog.pin_itp == prog_c.pin_itp
    @test prog.pout_itp == prog_c.pout_itp

    # Substance from the load-function-> 1st test
    db_path = string(@__DIR__, "/data")
    db = "Database_test.csv"
    solutes = ["C10", "C11"]
    init_t = zeros(length(solutes))
    init_τ = zeros(length(solutes))
    sub = GasChromatographySimulator.load_solute_database(db_path, db, sys.sp, sys.gas, solutes, init_t, init_τ)
    @test sub[1].name == "C10"
    @test sub[2].Tchar == 124.746 + 273.15

    # parameters
    par = GasChromatographySimulator.Parameters(sys, prog, sub, opt)
    @test par.sys.L == L
    @test par.prog.gf(0.0) == par.prog.gf(L) + par.prog.a_gf[:,1]
    @test par.prog.T_itp(0.0, 0.0) == temp_steps[1] + 273.15
    @test par.prog.T_itp(L, sum(time_steps)) == temp_steps[end] - ΔT_steps[end] + 273.15
    @test par.sub[1].Dag == GasChromatographySimulator.diffusivity(142.28, 10, 22, 0, 0, 0, "He")
    @test par.sub[2].Dag == GasChromatographySimulator.diffusivity(156.31, 11, 24, 0, 0, 0, "He")
end

@testset "solving check" begin
    # test simulation with different settings:
    # - with/without gradient
    # - Tcontrol = "inlet"/"outlet" (with gradient)
    # - odesys = true/false
    # - alg = OwrenZen3/OwrenZen5
    # same settings for:
    # - the System
    # - the Substances
    # - temperature and pressure program
    a_d = [d]
    a_df = [df]
    d_fun(x) = GasChromatographySimulator.gradient(x, a_d)
    df_fun(x) = GasChromatographySimulator.gradient(x, a_df)
    sys = GasChromatographySimulator.System(L, d_fun, a_d, df_fun, a_df, sp, gas)

    
    solutes = ["C10", "C11"]
    init_t = zeros(length(solutes))
    init_τ = zeros(length(solutes))
    sub = GasChromatographySimulator.load_solute_database(db_path, db, sys.sp, sys.gas, solutes, init_t, init_τ)

    opt = [ GasChromatographySimulator.Options(alg=OwrenZen3()),
            GasChromatographySimulator.Options(OwrenZen5(), 1e-6, 1e-3, "inlet", true),
            GasChromatographySimulator.Options(OwrenZen5(), 1e-6, 1e-3, "outlet", true),
            GasChromatographySimulator.Options(OwrenZen5(), 1e-6, 1e-3, "outlet", false),
            GasChromatographySimulator.Options(OwrenZen5(), 1e-6, 1e-3, "inlet", false),
            GasChromatographySimulator.Options(OwrenZen5(), 1e-8, 1e-5, "inlet", true)
            ]

    par_g = Array{GasChromatographySimulator.Parameters}(undef, length(opt))
    par_ng = Array{GasChromatographySimulator.Parameters}(undef, length(opt))
    results_g = Array{Any}(undef, length(opt))
    results_ng = Array{Any}(undef, length(opt))
    for i=1:length(opt)
        prog_g = GasChromatographySimulator.Program(time_steps, temp_steps, pin_steps, pout_steps, ΔT_steps, x₀_steps, L₀_steps, α_steps, opt[i].Tcontrol, sys.L)
        prog_ng = GasChromatographySimulator.Program(time_steps, temp_steps, pin_steps, pout_steps, [zeros(length(time_steps)) x₀_steps L₀_steps α_steps], opt[i].Tcontrol, sys.L)
        par_g[i] = GasChromatographySimulator.Parameters(sys, prog_g, sub, opt[i])
        par_ng[i] = GasChromatographySimulator.Parameters(sys, prog_ng, sub, opt[i])
        # simulate()
        results_g[i] = GasChromatographySimulator.simulate(par_g[i]) 
        results_ng[i] = GasChromatographySimulator.simulate(par_ng[i]) 
    end

    @test length(results_g[3]) == 2
    @test length(results_g[4]) == 3

    @test length(results_ng[3]) == 2
    @test length(results_ng[4]) == 3

    @test isapprox(results_g[1][1].tR[1], 123.184, atol=1e-3)
    @test isapprox(results_g[2][1].tR[1], 123.186, atol=1e-3)
    @test isapprox(results_g[3][1].tR[1], 51.4549, atol=1e-4)
    @test isapprox(results_g[4][1].tR[1], 51.4564, atol=1e-4)
    @test isapprox(results_g[5][1].tR[1], 123.208, atol=1e-3)
    @test isapprox(results_g[6][1].tR[1], 123.173, atol=1e-3)

    @test isapprox(results_g[1][1].τR[2], 0.548163, atol=1e-5)
    @test isapprox(results_g[2][1].τR[2], 0.547683, atol=1e-5)
    @test isapprox(results_g[3][1].τR[2], 0.525259, atol=1e-5)
    @test isapprox(results_g[4][1].τR[2], 0.524873, atol=1e-5)
    @test isapprox(results_g[5][1].τR[2], 0.547757, atol=1e-5)
    @test isapprox(results_g[6][1].τR[2], 0.547965, atol=1e-5)

    @test isapprox(results_g[1][1].Res[1], 18.5825, atol=1e-3)
    @test isapprox(results_g[2][1].Res[1], 18.5838, atol=1e-3)
    @test isapprox(results_g[3][1].Res[1], 19.8609, atol=1e-3)
    @test isapprox(results_g[4][1].Res[1], 19.8635, atol=1e-3)
    @test isapprox(results_g[5][1].Res[1], 18.5988, atol=1e-3)
    @test isapprox(results_g[6][1].Res[1], 18.5814, atol=1e-3)

    @test isapprox(results_ng[1][1].tR[1], 87.3973, atol=1e-3)
    @test isapprox(results_ng[2][1].tR[1], 87.4008, atol=1e-3)
    @test isapprox(results_ng[3][1].tR[1], 87.4008, atol=1e-3)
    @test isapprox(results_ng[4][1].tR[1], 87.4003, atol=1e-3)
    @test isapprox(results_ng[5][1].tR[1], 87.4003, atol=1e-3)
    @test isapprox(results_ng[6][1].tR[1], 87.3916, atol=1e-3)

    @test isapprox(results_ng[1][1].τR[2], 0.590034, atol=1e-5)
    @test isapprox(results_ng[2][1].τR[2], 0.590276, atol=1e-5)
    @test isapprox(results_ng[3][1].τR[2], 0.590276, atol=1e-5)
    @test isapprox(results_ng[4][1].τR[2], 0.595850, atol=1e-5)
    @test isapprox(results_ng[5][1].τR[2], 0.595850, atol=1e-5)
    @test isapprox(results_ng[6][1].τR[2], 0.590284, atol=1e-5)

    @test isapprox(results_ng[1][1].Res[1], 17.6590, atol=1e-3)
    @test isapprox(results_ng[2][1].Res[1], 17.6196, atol=1e-3)
    @test isapprox(results_ng[3][1].Res[1], 17.6196, atol=1e-3)
    @test isapprox(results_ng[4][1].Res[1], 17.5605, atol=1e-3)
    @test isapprox(results_ng[5][1].Res[1], 17.5605, atol=1e-3)
    @test isapprox(results_ng[6][1].Res[1], 17.6485, atol=1e-3)

    # sol_extraction()
    df_sol = GasChromatographySimulator.sol_extraction(results_g[1][2], par_g[1])
    @test df_sol.t[1][end] ≈ results_g[1][1].tR[1]
    @test df_sol.τ²[2][end] ≈ results_g[1][1].τR[2]^2

    df_sol = GasChromatographySimulator.sol_extraction(results_g[4][2], results_g[4][3], par_g[4])
    @test df_sol.t[1][end] ≈ results_g[4][1].tR[1]
    @test df_sol.τ²[2][end] ≈ results_g[4][1].τR[2]^2
end

@testset "misc checks" begin
    @test "C10" in GasChromatographySimulator.all_solutes(sp, DataFrame(CSV.File(string(db_path,"/",db), header=1, silencewarnings=true)))
    @test "C11" in GasChromatographySimulator.all_solutes(sp, DataFrame(CSV.File(string(db_path,"/",db), header=1, silencewarnings=true)))

    a_d = [d]
    a_df = [df]
    d_fun(x) = GasChromatographySimulator.gradient(x, a_d)
    df_fun(x) = GasChromatographySimulator.gradient(x, a_df)
    sys = GasChromatographySimulator.System(L, d_fun, a_d, df_fun, a_df, sp, gas)
    a_gf = [ΔT_steps x₀_steps L₀_steps α_steps]
    gf(x) = GasChromatographySimulator.gradient(x, a_gf)
    T_itp = GasChromatographySimulator.temperature_interpolation(time_steps, temp_steps, gf, sys.L)
    pin_itp = GasChromatographySimulator.pressure_interpolation(time_steps, pin_steps)
    pout_itp = GasChromatographySimulator.pressure_interpolation(time_steps, pout_steps)
    # without gradient
    prog_ng = GasChromatographySimulator.Program(time_steps, temp_steps, pin_steps, pout_steps, sys.L)

    # random time in the temperature program
    t = rand()*sum(time_steps)
    # temperature at random time at column inlet
    T_test = prog_ng.T_itp(0.0, t)

    η_T = GasChromatographySimulator.viscosity(T_test, "He")

    η_t = GasChromatographySimulator.viscosity(0.0, t, prog_ng.T_itp, sys.gas)

    @test η_T ≈ η_t

    tM_T = GasChromatographySimulator.holdup_time(T_test, prog_ng.pin_itp(t), prog_ng.pout_itp(t), sys.L, sys.a_d[1], sys.gas) # only defined for non-gradient case
    tM_t = GasChromatographySimulator.holdup_time(t, prog_ng.T_itp, prog_ng.pin_itp, prog_ng.pout_itp, sys.L, sys.d, sys.gas)

    @test tM_T ≈ tM_t

    F_T = GasChromatographySimulator.flow(T_test, prog_ng.pin_itp(t), prog_ng.pout_itp(t), sys.L, sys.a_d[1], sys.gas) # only defined for non-gradient case
    F_t = GasChromatographySimulator.flow(t, prog_ng.T_itp, prog_ng.pin_itp, prog_ng.pout_itp, sys.L, sys.d, sys.gas)

    @test F_T  ≈  F_t

end
#=
# solving
# for odesys:
sol_odesys = GasChromatographySimulator.solve_system_multithreads(par)
pl_odesys = GasChromatographySimulator.peaklist(sol_odesys, par)
# test: tR(C10) = 97.943, τR(C10) = 0.581
# test: tR(C11) = 140.416, τR(C11) = 0.540
@test round.(pl_odesys.tR, digits=3) == [97.943, 140.416]
@test round.(pl_odesys.τR, digits=3) == [0.581, 0.540]

opt_ng = GasChromatographySimulator.Options(ng=true)
par_ng = GasChromatographySimulator.Parameters(sys, prog_c_ng, sub, opt_ng)
sol_odesys_ng = GasChromatographySimulator.solve_system_multithreads(par_ng)
pl_odesys_ng = GasChromatographySimulator.peaklist(sol_odesys_ng, par_ng)
# non-gradient separation should give the 'same' results with both methods
# (ng=false and ng=true)
@test isapprox(pl_odesys_ng.tR, pl_odesys.tR; rtol=1e-6)
@test isapprox(pl_odesys_ng.τR, pl_odesys.τR; rtol=1e-4)

# for migration and peakvariance separatly:
sol_migr, sol_peak = GasChromatographySimulator.solve_multithreads(par)
pl_migr = GasChromatographySimulator.peaklist(sol_migr, sol_peak, par)
# test: tR(C10) = 97.944, τR(C10) = 0.578
# test: tR(C11) = 140.403, τR(C11) = 0.540
@test round.(pl_migr.tR, digits=3) == [97.944, 140.403]
@test round.(pl_migr.τR, digits=3) == [0.578, 0.540]

sol_migr_ng, sol_peak_ng = GasChromatographySimulator.solve_multithreads(par_ng)
pl_migr_ng = GasChromatographySimulator.peaklist(sol_migr_ng, sol_peak_ng, par_ng)
@test isapprox(pl_migr_ng.tR, pl_migr.tR; rtol=1e-5)
@test isapprox(pl_migr_ng.τR, pl_migr.τR; rtol=1e-4)

# test of holdup_time:
t = rand()*cumsum(time_steps)[end]
T_test = par.prog.T_itp(0.0, t)
η_T = GasChromatographySimulator.viscosity(T_test, "He")
η_t = GasChromatographySimulator.viscosity(0.0, t, par.prog.T_itp, par.sys.gas)

@test η_T == η_t

tM_T = GasChromatographySimulator.holdup_time(T_test, par.prog.pin_itp(t), par.prog.pout_itp(t), L, a_d[1], gas)
tM_t = GasChromatographySimulator.holdup_time(t, par.prog.T_itp, par.prog.pin_itp, par.prog.pout_itp, par.sys.L, par.sys.d, par.sys.gas; ng=false)
tM_t_ng = GasChromatographySimulator.holdup_time(t, par.prog.T_itp, par.prog.pin_itp, par.prog.pout_itp, par.sys.L, par.sys.d, par.sys.gas; ng=true)

@test tM_T ≈ tM_t
@test tM_t_ng ≈ tM_t

F_T = GasChromatographySimulator.flow(T_test, L, a_d[1], par.prog.pin_itp(t), par.prog.pout_itp(t), gas)
F_t = GasChromatographySimulator.flow(t, par.prog.T_itp, par.prog.pin_itp, par.prog.pout_itp, par.sys.L, par.sys.d, par.sys.gas; ng=false)
F_t_ng = GasChromatographySimulator.flow(t, par.prog.T_itp, par.prog.pin_itp, par.prog.pout_itp, par.sys.L, par.sys.d, par.sys.gas; ng=true)

# @test F_T  ≈  F_t
@test F_t_ng  ≈  F_t

# test with different tolerances 
opt_6_3 = GasChromatographySimulator.Options(OwrenZen5(), 1e-6, 1e-3, "inlet", true)
opt_6_4 = GasChromatographySimulator.Options(OwrenZen5(), 1e-6, 1e-4, "inlet", true)
opt_8_3 = GasChromatographySimulator.Options(OwrenZen5(), 1e-8, 1e-3, "inlet", true)
opt_8_4 = GasChromatographySimulator.Options(OwrenZen5(), 1e-8, 1e-4, "inlet", true)
# change somthing small, e.g. stretch the temperature program by a small nummber
# n = 1.000001 and compare the results for the same tolerances
# program with a gradient:
prog_0 = GasChromatographySimulator.constructor_Program(time_steps, temp_steps, pin_steps, pout_steps, 40.0.*ones(length(time_steps)), x₀_steps, L₀_steps, -3.0.*ones(length(time_steps)), opt.Tcontrol, sys.L)

par_6_3_0 = GasChromatographySimulator.Parameters(sys, prog_0, sub, opt_6_3)
par_6_4_0 = GasChromatographySimulator.Parameters(sys, prog_0, sub, opt_6_4)
par_8_3_0 = GasChromatographySimulator.Parameters(sys, prog_0, sub, opt_8_3)
par_8_4_0 = GasChromatographySimulator.Parameters(sys, prog_0, sub, opt_8_4)

pl_6_3_0 = GasChromatographySimulator.simulate(par_6_3_0)[1]
pl_6_4_0 = GasChromatographySimulator.simulate(par_6_4_0)[1]
pl_8_3_0 = GasChromatographySimulator.simulate(par_8_3_0)[1]
pl_8_4_0 = GasChromatographySimulator.simulate(par_8_4_0)[1]

prog_1 = GasChromatographySimulator.constructor_Program(1.000001.*time_steps, temp_steps, pin_steps, pout_steps, 40.0.*ones(length(time_steps)), x₀_steps, L₀_steps, -3.0.*ones(length(time_steps)), opt.Tcontrol, sys.L)
par_6_3_1 = GasChromatographySimulator.Parameters(sys, prog_1, sub, opt_6_3)
par_6_4_1 = GasChromatographySimulator.Parameters(sys, prog_1, sub, opt_6_4)
par_8_3_1 = GasChromatographySimulator.Parameters(sys, prog_1, sub, opt_8_3)
par_8_4_1 = GasChromatographySimulator.Parameters(sys, prog_1, sub, opt_8_4)

pl_6_3_1 = GasChromatographySimulator.simulate(par_6_3_1)[1]
pl_6_4_1 = GasChromatographySimulator.simulate(par_6_4_1)[1]
pl_8_3_1 = GasChromatographySimulator.simulate(par_8_3_1)[1]
pl_8_4_1 = GasChromatographySimulator.simulate(par_8_4_1)[1]

pl_6_3_0.tR .- pl_6_3_1.tR
pl_6_4_0.tR .- pl_6_4_1.tR
pl_8_3_0.tR .- pl_8_3_1.tR
pl_8_4_0.tR .- pl_8_4_1.tR

pl_6_3_0.tR .- pl_6_3_1.tR == pl_6_4_0.tR .- pl_6_4_1.tR
pl_6_3_0.tR .- pl_6_3_1.tR == pl_8_3_0.tR .- pl_8_3_1.tR
pl_6_3_0.tR .- pl_6_3_1.tR == pl_8_4_0.tR .- pl_8_4_1.tR
# how to interprete theses resultd
=#

println("Test run successful.")