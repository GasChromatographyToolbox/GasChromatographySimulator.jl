using Pkg
Pkg.activate("/Users/janleppert/Documents/GitHub/GasChromatographySimulator")
using GasChromatographySimulator
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
#@testset "parameters check" begin
    # default Options
    opt = GasChromatographySimulator.Options(OwrenZen5(), 1e-6, 1e-3, "inlet", true)
    GasChromatographySimulator.Options() == opt

    # System
    a_d = [d]
    a_df = [df]
    d_fun(x) = GasChromatographySimulator.gradient(x, a_d)
    df_fun(x) = GasChromatographySimulator.gradient(x, a_df)
    sys = GasChromatographySimulator.System(10.0, d_fun, a_d, df_fun, a_df, sp, gas)
    sys_c = GasChromatographySimulator.constructor_System(L, a_d[1], a_df[1], sp, gas)
    sys.d(sys.L) == sys_c.d(sys_c.L)

    # Program 
    a_gf = [ΔT_steps x₀_steps L₀_steps α_steps]
    gf(x) = GasChromatographySimulator.gradient(x, a_gf)
    T_itp = GasChromatographySimulator.temperature_interpolation(time_steps, temp_steps, gf, sys.L)
    pin_itp = GasChromatographySimulator.pressure_interpolation(time_steps, pin_steps)
    pout_itp = GasChromatographySimulator.pressure_interpolation(time_steps, pout_steps)
    prog = GasChromatographySimulator.Program(time_steps, temp_steps, pin_steps, pout_steps, gf, a_gf, T_itp, pin_itp, pout_itp)
    prog_c = GasChromatographySimulator.constructor_Program(time_steps, temp_steps, pin_steps, pout_steps, ΔT_steps, x₀_steps, L₀_steps, α_steps, opt.Tcontrol, sys.L)
    prog_c_ng = GasChromatographySimulator.constructor_Program(time_steps, temp_steps, pin_steps, pout_steps, sys.L) # without gradient
    prog.T_itp == prog_c.T_itp
    prog_c.T_itp != prog_c_ng.T_itp
    prog.pin_itp == prog_c.pin_itp
    prog.pout_itp == prog_c.pout_itp

    # Substance from the load-function-> 1st test
    db_path = string(@__DIR__, "/data")
    db = "Database_test.csv"
    solutes = ["C10", "C11"]
    init_t = zeros(length(solutes))
    init_τ = zeros(length(solutes))
    sub = GasChromatographySimulator.load_solute_database(db_path, db, sys.sp, sys.gas, solutes, init_t, init_τ)
    sub[1].name == "C10"
    sub[2].Tchar == 124.746 + 273.15

    # parameters
    par = GasChromatographySimulator.Parameters(sys, prog, sub, opt)
    par.sys.L == L
    par.prog.gf(0.0) == par.prog.gf(L) + par.prog.a_gf[:,1]
    par.prog.T_itp(0.0, 0.0) == temp_steps[1] + 273.15
    par.prog.T_itp(L, sum(time_steps)) == temp_steps[end] - ΔT_steps[end] + 273.15
    par.sub[1].Dag == GasChromatographySimulator.diffusivity(142.28, 10, 22, 0, 0, 0, "He")
    par.sub[2].Dag == GasChromatographySimulator.diffusivity(156.31, 11, 24, 0, 0, 0, "He")
#end

#@testset "solving and results check" begin
    # test simulation with different settings:
    # - with/without gradient
    # - Tcontrol = "inlet"/"outlet" (with gradient)
    # - odesys = true/false
    # - alg = OwrenZen3/OwrenZen5
    # - outlet pressure = 101.3kPa/0.0kPa
    # same settings for:
    # - the System
    # - the Substances
    # - temperature and pressure program
    a_d = [d]
    a_df = [df]
    d_fun(x) = GasChromatographySimulator.gradient(x, a_d)
    df_fun(x) = GasChromatographySimulator.gradient(x, a_df)
    sys = GasChromatographySimulator.System(10.0, d_fun, a_d, df_fun, a_df, sp, gas)

    db_path = string(@__DIR__, "/data")
    db = "Database_test.csv"
    solutes = ["C10", "C11"]
    init_t = zeros(length(solutes))
    init_τ = zeros(length(solutes))
    sub = GasChromatographySimulator.load_solute_database(db_path, db, sys.sp, sys.gas, solutes, init_t, init_τ)

    opt = [ GasChromatographySimulator.Options(OwrenZen3(), 1e-6, 1e-3, "inlet", true),
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
        prog_ng = GasChromatographySimulator.Program(time_steps, temp_steps, pin_steps, pout_steps, zeros(length(time_steps)), x₀_steps, L₀_steps, α_steps, opt[i].Tcontrol, sys.L)
        par_g[i] = GasChromatographySimulator.Parameters(sys, prog_g, sub, opt[i])
        par_ng[i] = GasChromatographySimulator.Parameters(sys, prog_ng, sub, opt[i])
        # simulate()
        results_g[i] = GasChromatographySimulator.simulate(par_g[i]) 
        results_ng[i] = GasChromatographySimulator.simulate(par_ng[i]) 
    end

    results_g[1][1].tR[1] == 123.18424848999621
    results_g[2][1].tR[1] == 123.18624737819937
    results_g[3][1].tR[1] == 51.454863869253934
    results_g[4][1].tR[1] == 51.45638645568056
    results_g[5][1].tR[1] == 123.20828349634378
    results_g[6][1].tR[1] == 123.17339663387685

    results_g[1][1].τR[2] == 0.5481628763675231
    results_g[2][1].τR[2] == 0.5476833396406813
    results_g[3][1].τR[2] == 0.5252588001145008
    results_g[4][1].τR[2] == 0.5248726201963742
    results_g[5][1].τR[2] == 0.5477573346176421
    results_g[6][1].τR[2] == 0.5479648788481386

    results_g[1][1].Res[1] == 18.58246263731392
    results_g[2][1].Res[1] == 18.583767746802923
    results_g[3][1].Res[1] == 19.860914280255994
    results_g[4][1].Res[1] == 19.863463267876234
    results_g[5][1].Res[1] == 18.598827557069868
    results_g[6][1].Res[1] == 18.581356452824906

    results_ng[1][1].tR[1] == 87.39733183107204
    results_ng[2][1].tR[1] == 87.40081993839854
    results_ng[3][1].tR[1] == 87.40081993839854
    results_ng[4][1].tR[1] == 87.40030480193386
    results_ng[5][1].tR[1] == 87.40030480193386
    results_ng[6][1].tR[1] == 87.39160123672703

    results_ng[1][1].τR[2] == 0.5900336247965692
    results_ng[2][1].τR[2] == 0.5902763609009182
    results_ng[3][1].τR[2] == 0.5902763609009182
    results_ng[4][1].τR[2] == 0.5958501556434311
    results_ng[5][1].τR[2] == 0.5958501556434311
    results_ng[6][1].τR[2] == 0.5902842414446662

    results_ng[1][1].Res[1] == 17.65902199872825
    results_ng[2][1].Res[1] == 17.619585565883256
    results_ng[3][1].Res[1] == 17.619585565883256
    results_ng[4][1].Res[1] == 17.560464952423732
    results_ng[5][1].Res[1] == 17.560464952423732
    results_ng[6][1].Res[1] == 17.648531928432124

    df_sol = GasChromatographySimulator.sol_extraction(results_g[1][2], par_g[1])
    df_sol.t[1][end] == results_g[1][1].tR[1]
    df_sol.τ²[2][end] ≈ results_g[1][1].τR[2]^2

    df_sol = GasChromatographySimulator.sol_extraction(results_g[4][2], results_g[4][3], par_g[4])
    df_sol.t[1][end] == results_g[4][1].tR[1]
    df_sol.τ²[2][end] ≈ results_g[4][1].τR[2]^2
#end

# misc testset

"C10" in GasChromatographySimulator.all_solutes(sp, DataFrame(CSV.File(string(db_path,"/",db), header=1, silencewarnings=true)))

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
    # with gradient
    prog = GasChromatographySimulator.Program(time_steps, temp_steps, pin_steps, pout_steps, gf, a_gf, T_itp, pin_itp, pout_itp)
    # without gradient
    prog_ng = GasChromatographySimulator.Program(time_steps, temp_steps, pin_steps, pout_steps, sys.L)

    # random time in the temperature program
    t = rand()*sum(time_steps)
    # temperature at random time at column inlet
    T_test = prog_ng.T_itp(0.0, t)

    η_T = GasChromatographySimulator.viscosity(T_test, "He")

    η_t = GasChromatographySimulator.viscosity(0.0, t, prog_ng.T_itp, sys.gas)

    η_T ≈ η_t

    tM_T = GasChromatographySimulator.holdup_time(T_test, prog_ng.pin_itp(t), prog_ng.pout_itp(t), sys.L, sys.a_d[1], sys.gas) # only defined for non-gradient case
    tM_t = GasChromatographySimulator.holdup_time(t, prog_ng.T_itp, prog_ng.pin_itp, prog_ng.pout_itp, sys.L, sys.d, sys.gas)

    tM_T ≈ tM_t

    F_T = GasChromatographySimulator.flow(T_test, prog_ng.pin_itp(t), prog_ng.pout_itp(t), sys.L, sys.a_d[1], sys.gas) # only defined for non-gradient case
    F_t = GasChromatographySimulator.flow(t, prog_ng.T_itp, prog_ng.pin_itp, prog_ng.pout_itp, sys.L, sys.d, sys.gas)

    F_T  ≈  F_t