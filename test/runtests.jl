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
            GasChromatographySimulator.Options(OwrenZen5(), 1e-6, 1e-3, "inlet", false)#,
            #GasChromatographySimulator.Options(OwrenZen5(), 1e-6, 1e-3, "inlet", true; ng=true)
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
    #@test isapprox(results_g[6][1].tR[1], 123.173, atol=1e-3)

    @test isapprox(results_g[1][1].τR[2], 0.548163, atol=1e-5)
    @test isapprox(results_g[2][1].τR[2], 0.547683, atol=1e-5)
    @test isapprox(results_g[3][1].τR[2], 0.525259, atol=1e-5)
    @test isapprox(results_g[4][1].τR[2], 0.524873, atol=1e-4)
    @test isapprox(results_g[5][1].τR[2], 0.547757, atol=1e-5)
    #@test isapprox(results_g[6][1].τR[2], 0.547965, atol=1e-5)

    @test isapprox(results_g[1][1].Res[1], 18.5825, atol=1e-3)
    @test isapprox(results_g[2][1].Res[1], 18.5838, atol=1e-3)
    @test isapprox(results_g[3][1].Res[1], 19.8609, atol=1e-3)
    @test isapprox(results_g[4][1].Res[1], 19.8635, atol=1e-3)
    @test isapprox(results_g[5][1].Res[1], 18.5988, atol=1e-3)
    #@test isapprox(results_g[6][1].Res[1], 18.5814, atol=1e-3)

    @test isapprox(results_ng[1][1].tR[1], 87.3973, atol=1e-3)
    @test isapprox(results_ng[2][1].tR[1], 87.4008, atol=1e-3)
    @test isapprox(results_ng[3][1].tR[1], 87.4008, atol=1e-3)
    @test isapprox(results_ng[4][1].tR[1], 87.4003, atol=1e-3)
    @test isapprox(results_ng[5][1].tR[1], 87.4003, atol=1e-3)
    #@test isapprox(results_ng[6][1].tR[1], 87.3916, atol=1e-3)

    @test isapprox(results_ng[1][1].τR[2], 0.590034, atol=1e-5)
    @test isapprox(results_ng[2][1].τR[2], 0.590276, atol=1e-5)
    @test isapprox(results_ng[3][1].τR[2], 0.590276, atol=1e-5)
    @test isapprox(results_ng[4][1].τR[2], 0.595850, atol=1e-5)
    @test isapprox(results_ng[5][1].τR[2], 0.595850, atol=1e-5)
    #@test isapprox(results_ng[6][1].τR[2], 0.590284, atol=1e-5)

    @test isapprox(results_ng[1][1].Res[1], 17.6590, atol=1e-3)
    @test isapprox(results_ng[2][1].Res[1], 17.6196, atol=1e-3)
    @test isapprox(results_ng[3][1].Res[1], 17.6196, atol=1e-3)
    @test isapprox(results_ng[4][1].Res[1], 17.5605, atol=1e-3)
    @test isapprox(results_ng[5][1].Res[1], 17.5605, atol=1e-3)
    #@test isapprox(results_ng[6][1].Res[1], 17.6485, atol=1e-3)

    # sol_extraction()
    df_sol = GasChromatographySimulator.sol_extraction(results_g[1][2], par_g[1])
    @test df_sol.t[1][end] ≈ results_g[1][1].tR[1]
    @test df_sol.τ²[2][end] ≈ results_g[1][1].τR[2]^2

    df_sol = GasChromatographySimulator.sol_extraction(results_g[4][2], results_g[4][3], par_g[4])
    @test df_sol.t[1][end] ≈ results_g[4][1].tR[1]
    @test df_sol.τ²[2][end] ≈ results_g[4][1].τR[2]^2

    # ng = true
    opt_ng = GasChromatographySimulator.Options(ng=true)
    prog_o_ng = GasChromatographySimulator.Program(time_steps, temp_steps, pin_steps, pout_steps, [zeros(length(time_steps)) x₀_steps L₀_steps α_steps], opt_ng.Tcontrol, sys.L)
    par_o_ng = GasChromatographySimulator.Parameters(sys, prog_o_ng, sub, opt_ng)
    results_o_ng = GasChromatographySimulator.simulate(par_o_ng)
    @test isapprox(results_o_ng[1].tR[1], 87.401, atol=1e-3)
    @test isapprox(results_o_ng[1].τR[2], 0.59028, atol=1e-5)

    # vis = "HP"
    opt_vis = GasChromatographySimulator.Options(vis="HP")
    prog_vis = GasChromatographySimulator.Program(time_steps, temp_steps, pin_steps, pout_steps, [zeros(length(time_steps)) x₀_steps L₀_steps α_steps], opt_vis.Tcontrol, sys.L)
    par_vis = GasChromatographySimulator.Parameters(sys, prog_vis, sub, opt_vis)
    results_vis = GasChromatographySimulator.simulate(par_vis)
    isapprox(results_vis[1].tR[1], 87.078, atol=1e-3)
    isapprox(results_vis[1].τR[2], 0.59294, atol=1e-5)
end



@testset "plots check" begin

    a_d = [d]
    a_df = [df]
    d_fun(x) = GasChromatographySimulator.gradient(x, a_d)
    df_fun(x) = GasChromatographySimulator.gradient(x, a_df)
    sys = GasChromatographySimulator.System(L, d_fun, a_d, df_fun, a_df, sp, gas)

    solutes = ["C10", "C11"]
    init_t = zeros(length(solutes))
    init_τ = zeros(length(solutes))
    sub = GasChromatographySimulator.load_solute_database(db_path, db, sys.sp, sys.gas, solutes, init_t, init_τ)
    
    opt_ng = GasChromatographySimulator.Options(ng=true)
    prog_o_ng = GasChromatographySimulator.Program(time_steps, temp_steps, pin_steps, pout_steps, [zeros(length(time_steps)) x₀_steps L₀_steps α_steps], opt_ng.Tcontrol, sys.L)
    par_o_ng = GasChromatographySimulator.Parameters(sys, prog_o_ng, sub, opt_ng)
    results_o_ng = GasChromatographySimulator.simulate(par_o_ng)

    # separat testset
    # plot functions
    t_start = 0.0
    t_end = 200.0
    p_chrom, t, chrom = GasChromatographySimulator.plot_chromatogram(results_o_ng[1], (t_start, t_end))
    @test t[1] == t_start && t[end] == t_end
    tR = results_o_ng[1].tR[1]
    τR = results_o_ng[1].τR[1]
    t0 = round(tR; digits=1)
    @test chrom[findfirst(t.==t0)] == 1/sqrt(2*π*τR^2)*exp(-(t0-tR)^2/(2*τR^2))
    t2, chrom2 = GasChromatographySimulator.plot_chromatogram!(p_chrom, results_o_ng[1], (t_start, t_end); mirror=true)
    @test t == t2
    @test chrom == -chrom2
    #@test p_chrom[1][:xaxis][:optimized_ticks][1][1] = t_start 
    #@test p_chrom[1][:xaxis][:optimized_ticks][1][end] = t_end 
    @test p_chrom[1][1][:y] == - p_chrom[1][2][:y]
    @test p_chrom[1][1][:x] == p_chrom[1][2][:x]

    p_flow = GasChromatographySimulator.plot_flow(par_o_ng)
    #@test p_flow[1][:xaxis][:optimized_ticks][1][1] = 0.0 
    #@test p_flow[1][:xaxis][:optimized_ticks][1][end] = 800.0
    @test isapprox(p_flow[1][1][:y][1], 11.9797, atol=1e-4)

    p_pres = GasChromatographySimulator.plot_pressure(par_o_ng.prog)
    #@test p_pres[1][:xaxis][:optimized_ticks][1][1] = 0.0 
    #@test p_pres[1][:xaxis][:optimized_ticks][1][end] = 800.0
    @test p_pres[1][1][:y][1] == par_o_ng.prog.pin_steps[1]
    @test p_pres[1][2][:y][end] == par_o_ng.prog.pout_steps[end]

    p_temp = GasChromatographySimulator.plot_temperature(par_o_ng)
    @test p_temp[1][1][:y][1] == par_o_ng.prog.temp_steps[1]
    @test p_temp[1][1][:y][end] == par_o_ng.prog.temp_steps[end]
    @test p_temp[1][2][:x][end] == sum(par_o_ng.prog.time_steps)

    p_temp = GasChromatographySimulator.plot_temperature(par_o_ng; selector="T(x)")
    @test p_temp[1][1][:y][1] == par_o_ng.prog.temp_steps[1]
    @test p_temp[1][end][:y][end] == par_o_ng.prog.temp_steps[end]
    @test p_temp[1][2][:x][end] == par_o_ng.sys.L

    p_temp = GasChromatographySimulator.plot_temperature(par_o_ng; selector="T(x,t)")
    @test p_temp[1][1][:y][end] == sum(par_o_ng.prog.time_steps)
    @test p_temp[1][1][:x][end] == par_o_ng.sys.L
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

println("Test run successful.")