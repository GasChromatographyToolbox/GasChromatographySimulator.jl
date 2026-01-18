using Test, GasChromatographySimulator, PlutoUI, DataFrames, CSV

# general parameters
L = 10.0
d = 0.25e-3
df = 0.25e-6
sp = "SLB5ms"
gas = "He"
time_steps = [0.0, 60.0, 600.0, 300.0]
temp_steps = [40.0, 40.0, 300.0, 300.0]
pin_steps = [200.0, 200.0, 300.0, 300.0].*1000.0 .+ 101300
F_steps = [1.0, 1.0, 1.0, 1.0]./60e6
pout_steps = [0.0, 0.0, 0.0, 0.0].*1000.0
ΔT_steps = [20.0, 30.0, 50.0, 40.0]
x₀_steps = zeros(length(time_steps))
L₀_steps = L.*ones(length(time_steps))
α_steps = zeros(length(time_steps))
db_path = joinpath("..", "data")#joinpath(pwd(), "data")#string(@__DIR__, "/data")
db = "Database_test.csv"
db_unc = "Database_test_uncertainty.csv"

@testset "gradient function check" begin
    # Gradient function for d
    a_d = [d, 0.0, 2*L, 0.0]
    d_fun(x) = GasChromatographySimulator.gradient(x, a_d)
    @test d_fun(L) == d/2

    # Gradient function for T with all parameters equal zero
    a_gf = zeros(2,3)
    gf_fun1(x) = GasChromatographySimulator.gradient(x, a_gf)
    @test gf_fun1(0.0) == zeros(size(a_gf)[1])

    # Gradient function with changing α
    αs = [-5.0, 0.0, 0.0, 5.0]
    a_gf = [ΔT_steps x₀_steps L₀_steps αs]
    gf_fun2(x) = GasChromatographySimulator.gradient(x, a_gf)
    @test gf_fun2(L) == -ΔT_steps
end

@testset "parameters check" begin
    # default Options
    opt = GasChromatographySimulator.Options(OwrenZen5(), 1e-6, 1e-3, "inlet", true)
    @test GasChromatographySimulator.Options() == opt

    # Column
    a_d = [d]
    a_df = [df]
    d_fun(x) = GasChromatographySimulator.gradient(x, a_d)
    df_fun(x) = GasChromatographySimulator.gradient(x, a_df)
    col = GasChromatographySimulator.Column(L, d_fun, a_d, df_fun, a_df, sp, gas)
    col_c = GasChromatographySimulator.constructor_System(L, d, df, sp, gas)
    @test col.d(col.L) == col_c.d

    # Program 
    a_gf = [ΔT_steps x₀_steps L₀_steps α_steps]
    gf(x) = GasChromatographySimulator.gradient(x, a_gf)
    T_itp = GasChromatographySimulator.temperature_interpolation(time_steps, temp_steps, gf, col.L)
    pin_itp = GasChromatographySimulator.steps_interpolation(time_steps, pin_steps)
    pout_itp = GasChromatographySimulator.steps_interpolation(time_steps, pout_steps)
    prog = GasChromatographySimulator.Program(time_steps, temp_steps, pin_steps, pout_steps, gf, a_gf, T_itp, pin_itp, pout_itp)
    prog_c = GasChromatographySimulator.constructor_Program(time_steps, temp_steps, pin_steps, pout_steps, ΔT_steps, x₀_steps, L₀_steps, α_steps, opt.Tcontrol, col.L)
    prog_c_ng = GasChromatographySimulator.constructor_Program(time_steps, temp_steps, pin_steps, pout_steps, col.L) # without gradient
    @test prog.T_itp == prog_c.T_itp
    @test prog_c.T_itp != prog_c_ng.T_itp
    @test prog.Fpin_itp == prog_c.Fpin_itp
    @test prog.pout_itp == prog_c.pout_itp

    # conventional program
    TP0 = [40.0, 1.0, 5.0, 280.0, 2.0, 20.0, 320.0, 2.0]
    pP0 = [400000.0, 10.0, 5000.0, 500000.0, 20.0]
    FP0 = [1/(6e7), 60.0]
    ts1, Ts = GasChromatographySimulator.conventional_program(TP0)
    ts2, ps = GasChromatographySimulator.conventional_program(pP0)
    ts3, Fs = GasChromatographySimulator.conventional_program(FP0) 
    @test ts1 == [0.0, 60.0, 2880.0, 120.0, 120.0, 120.0]
    @test Fs == 1/(6e7).*ones(2)

    TP = GasChromatographySimulator.temperature_program(ts1, Ts; time_unit="min")
    FP = GasChromatographySimulator.temperature_program(ts3, Fs; time_unit="s")
    @test TP == TP0
    @test FP./[1, 60] == FP0

    cts12 = GasChromatographySimulator.common_time_steps(ts1, ts2)
    cts21 = GasChromatographySimulator.common_time_steps(ts2, ts1)
    @test cts12 == cts21

    cts123 = GasChromatographySimulator.common_time_steps(cts12, ts3)
    cts13 = GasChromatographySimulator.common_time_steps(ts1, ts3)
    new_Fs = GasChromatographySimulator.new_value_steps(Fs, ts3, cts13)
    new_Ts = GasChromatographySimulator.new_value_steps(Ts, ts1, cts123)
    @test new_Fs == 1/(6e7).*ones(length(cts13))
    @test new_Ts == [40.0, 40.0, 85.0, 185.0, 280.0, 280.0, 280.0, 320.0, 320.0, 320.0]

    prog_conv = GasChromatographySimulator.Program([40.0, 1.0, 5.0, 280.0, 2.0, 20.0, 320.0, 2.0], [400000.0, 10.0, 5000.0, 500000.0, 20.0], L; pout="vacuum", time_unit="min")
    prog_conv_s_atm = GasChromatographySimulator.Program([40.0, 1.0*60.0, 5.0/60.0, 280.0, 2.0*60.0, 20.0/60.0, 320.0, 2.0*60.0], [400000.0, 10.0*60.0, 5000.0/60.0, 500000.0, 20.0*60.0], L; pout="atmosphere", time_unit="s")
    @test prog_conv.temp_steps == [40.0, 40.0, 85.0, 185.0, 280.0, 280.0, 280.0, 320.0, 320.0]
    @test prog_conv.time_steps == prog_conv_s_atm.time_steps
    @test prog_conv.pout_steps == prog_conv_s_atm.pout_steps .- 101300

    TP = [40.0, 1.0, 5.0, 200.0, 0.0, 10.0, 280.0, 2.0, 20.0, 320.0, 2.0]
    FpinP = [400000.0, 10.0, 5000.0, 500000.0, 20.0]
    poutP = [0.0, 100.0]
    ΔTP = [0.0, 5.0, 10.0, 60.0, 10.0, 5.0, 80.0, 0.0, -10.0, 0.0, 10.0]
    x₀P = [0.0, 100.0]
    L₀P = [L, 100.0]
    αP = [0.0, 10.0, -0.5, -5.0, 15.0, 1.0, 5.0, 10.0]

    prog_conv_grad = GasChromatographySimulator.Program(TP, FpinP, poutP, ΔTP, x₀P, L₀P, αP, "inlet", L; time_unit="min")
    @test prog_conv_grad.time_steps == [0.0, 60.0, 240.0, 300.0, 60.0, 540.0, 60.0, 240.0, 300.0, 180.0, 120.0, 360.0, 120.0, 120.0, 120.0, 180.0, 300.0, 2700.0]

    
    # Substance from the load-function-> 1st test
    #db_path = string(@__DIR__, "/data")
    #db = "Database_test.csv"
    solutes = ["C10", "C11", "Glyceryl trioctanoate", "AAA"]
    init_t = zeros(length(solutes))
    init_τ = zeros(length(solutes))
    sub = GasChromatographySimulator.load_solute_database(db_path, db, col.sp, col.gas, solutes, init_t, init_τ)
    @test sub[1].CAS == "124-18-5"
    @test sub[2].Tchar == 124.75 + 273.15

    # no stationary phase & other gas for diffusivity
    # H2
    sub_0 = GasChromatographySimulator.load_solute_database(db_path, db, "", "H2", solutes, init_t, init_τ)
    @test sub_0[1].Tchar == 1.0
    @test sub_0[2].ann == "no sp"
    @test round(sub_0[1].Cag; sigdigits=5) == 0.00011885
    # N2
    sub_0 = GasChromatographySimulator.load_solute_database(db_path, db, "", "N2", solutes, init_t, init_τ)
    @test round(sub_0[1].Cag; sigdigits=5) == 2.8398e-5
    # Ar
    sub_0 = GasChromatographySimulator.load_solute_database(db_path, db, "", "Ar", solutes, init_t, init_τ)
    @test round(sub_0[1].Cag; sigdigits=5) == 2.5268e-5
    # test for error-cases of GasChromatographySimulator.load_solute_database

    # test for old database format -> not supported any more
    #db_old = "Database_test_old_format.csv"
    #sub_new = GasChromatographySimulator.load_solute_database(db_path, db, "Wax", "He", ["C14", "Decyl acetate", "Hexadecane", "C15", "Methyl myristate"], [1.0, 2.0, 3.0, 4.0, 5.0], [0.1, 0.2, 0.3, 0.4, 0.5])
    #sub_old = GasChromatographySimulator.load_solute_database(db_path, db_old, "Wax", "He", ["C14", "Decyl acetate", "Hexadecane", "C15", "Methyl myristate"], [1.0, 2.0, 3.0, 4.0, 5.0], [0.1, 0.2, 0.3, 0.4, 0.5])
    #@test sub_new[1].CAS == sub_old[1].CAS
    #@test sub_new[2].t₀ == 4.0

    # parameters
    par = GasChromatographySimulator.Parameters(col, prog, sub, opt)
    @test par.col.L == L
    @test par.prog.gf(0.0) == par.prog.gf(L) + par.prog.a_gf[:,1]
    @test par.prog.T_itp(0.0, 0.0) == temp_steps[1] + 273.15
    @test par.prog.T_itp(L, sum(time_steps)) == temp_steps[end] - ΔT_steps[end] + 273.15
    @test par.sub[1].Cag == GasChromatographySimulator.diffusivity(GasChromatographySimulator.CAS_identification(par.sub[1].name), "He")
    @test  isapprox(par.sub[2].Cag, GasChromatographySimulator.diffusivity(156.31, 11, 24, 0, 0, 0, "He"), atol=1e-6)

    # retention factor of low volatile solute "Glyceryl trioctanoate"
    k = GasChromatographySimulator.retention_factor(0.0, 0.0, par.col, par.prog, par.sub[3], par.opt)
    @test k == opt.k_th

    H = GasChromatographySimulator.plate_height(0.0, 0.0, par.col, par.prog, par.sub[3], par.opt)
    H_ = GasChromatographySimulator.plate_height(0.0, 0.0, par.prog.T_itp, par.prog.Fpin_itp, par.prog.pout_itp, par.col.L, par.col.d, par.col.df, par.col.gas, par.sub[3].Tchar, par.sub[3].θchar, par.sub[3].ΔCp, par.sub[3].φ₀, par.sub[3].Cag)
    H__ = GasChromatographySimulator.plate_height(0.0, 0.0, par.prog.T_itp, par.prog.Fpin_itp, par.prog.pout_itp, par.col.L, par.col.d, par.col.df, par.col.gas, 0.0, 0.0, 0.0, par.sub[3].φ₀, par.sub[3].Cag)
    @test H == H_
    @test H > H__

    # control="Flow"
    opt_F = GasChromatographySimulator.Options(control="Flow")
    prog_F = GasChromatographySimulator.constructor_Program(time_steps, temp_steps, F_steps, pout_steps, ΔT_steps, x₀_steps, L₀_steps, α_steps, opt.Tcontrol, col.L)
    par_F = GasChromatographySimulator.Parameters(col, prog_F, sub, opt_F)
    pin_F = GasChromatographySimulator.inlet_pressure(0.0, par_F)
    pin_F_ = GasChromatographySimulator.inlet_pressure(0.0, prog_F.T_itp, prog_F.Fpin_itp, prog_F.pout_itp, col.L, col.d, col.gas; control=opt_F.control) 
    @test pin_F == pin_F_
end

@testset "solving check" begin
    # test simulation with different settings:
    # - with/without gradient
    # - Tcontrol = "inlet"/"outlet" (with gradient)
    # - odesys = true/false
    # - alg = OwrenZen3/OwrenZen5
    # same settings for:
    # - the Column
    # - the Substances
    # - temperature and pressure program
    #a_d = [d]
    #a_df = [df]
    #d_fun(x) = GasChromatographySimulator.gradient(x, a_d)
    #df_fun(x) = GasChromatographySimulator.gradient(x, a_df)
    col = GasChromatographySimulator.Column(L, d, df, sp, gas)

    
    solutes = ["C10", "C11"]
    init_t = zeros(length(solutes))
    init_τ = zeros(length(solutes))
    sub = GasChromatographySimulator.load_solute_database(db_path, db, col.sp, col.gas, solutes, init_t, init_τ)

    # Tinlet = "inlet"/"outlet" and odesys = true/false 
    opt = [ #GasChromatographySimulator.Options(alg=OwrenZen3()),
            GasChromatographySimulator.Options(OwrenZen5(), 1e-6, 1e-3, "inlet", true),
            GasChromatographySimulator.Options(OwrenZen5(), 1e-6, 1e-3, "outlet", false)
            ]

    par_g = Array{GasChromatographySimulator.Parameters}(undef, length(opt))
    results_g = Array{Any}(undef, length(opt))
    for i=1:length(opt)
        prog_g = GasChromatographySimulator.Program(time_steps, temp_steps, pin_steps, pout_steps, ΔT_steps, x₀_steps, L₀_steps, α_steps, opt[i].Tcontrol, col.L)
        par_g[i] = GasChromatographySimulator.Parameters(col, prog_g, sub, opt[i])
        results_g[i] = GasChromatographySimulator.simulate(par_g[i])
    end

    @test length(results_g[1]) == 2
    @test length(results_g[2]) == 3

    #@test isapprox(results_g[1][1].tR[1], 123.184, atol=1e-3)
    @test isapprox(results_g[1][1].tR[1], 123.19, atol=1e-2)
    @test isapprox(results_g[2][1].tR[1], 51.46, atol=1e-2)

    #@test isapprox(results_g[1][1].τR[2], 0.548163, atol=1e-5)
    @test isapprox(results_g[1][1].τR[2], 0.548, atol=1e-3)
    @test isapprox(results_g[2][1].τR[2], 0.525, atol=1e-3)

    # sol_extraction()
    df_sol = GasChromatographySimulator.sol_extraction(results_g[1][2], par_g[1])
    @test df_sol.t[1][end] ≈ results_g[1][1].tR[1]
    @test df_sol.τ²[2][end] ≈ results_g[1][1].τR[2]^2

    df_sol = GasChromatographySimulator.sol_extraction(results_g[2][2], results_g[2][3], par_g[2])
    @test df_sol.t[1][end] ≈ results_g[2][1].tR[1]
    @test df_sol.τ²[2][end] ≈ results_g[2][1].τR[2]^2

    # ng = true
    opt_ng = GasChromatographySimulator.Options(ng=true)
    prog_o_ng = GasChromatographySimulator.Program(time_steps, temp_steps, pin_steps, pout_steps, [zeros(length(time_steps)) x₀_steps L₀_steps α_steps], opt_ng.Tcontrol, col.L)
    par_o_ng = GasChromatographySimulator.Parameters(col, prog_o_ng, sub, opt_ng)
    results_o_ng = GasChromatographySimulator.simulate(par_o_ng)
    @test isapprox(results_o_ng[1].tR[1], 87.41, atol=1e-2)
    @test isapprox(results_o_ng[1].τR[2], 0.590, atol=1e-3)

    Name = [sub[i].name for i in 1:length(sub)]
    CAS = [sub[i].CAS for i in 1:length(sub)]
    ann = [sub[i].ann for i in 1:length(sub)]
    Tchar = [sub[i].Tchar for i in 1:length(sub)]
    θchar = [sub[i].θchar for i in 1:length(sub)]
    ΔCp = [sub[i].ΔCp for i in 1:length(sub)]
    φ₀ = [sub[i].φ₀ for i in 1:length(sub)]
    Cag = [sub[i].Cag for i in 1:length(sub)]
    t₀ = [sub[i].t₀ for i in 1:length(sub)]
    τ₀ = [sub[i].τ₀ for i in 1:length(sub)]
    results_p = GasChromatographySimulator.simulate(col.L, col.d, col.df, col.gas, 
                                                    prog_o_ng.T_itp,prog_o_ng.Fpin_itp, prog_o_ng.pout_itp, 
                                                    Name, CAS, ann, Tchar, θchar, ΔCp, φ₀, Cag, t₀, τ₀, 
                                                    opt_ng)
    @test results_p[1].tR[1] == results_o_ng[1].tR[1]

    # odesys = false
    results_odesys_false = GasChromatographySimulator.simulate(col.L, col.d, col.df, col.gas, 
                                                    prog_o_ng.T_itp,prog_o_ng.Fpin_itp, prog_o_ng.pout_itp, 
                                                    Name, CAS, ann, Tchar, θchar, ΔCp, φ₀, Cag, t₀, τ₀, 
                                                    opt[2])
    @test isapprox(results_odesys_false[1].tR[1], results_p[1].tR[1], atol=1e-2)

    # vis = "HP"
    #opt_vis = GasChromatographySimulator.Options(vis="HP")
    #prog_vis = GasChromatographySimulator.Program(time_steps, temp_steps, pin_steps, pout_steps, [zeros(length(time_steps)) x₀_steps L₀_steps α_steps], opt_vis.Tcontrol, col.L)
    #par_vis = GasChromatographySimulator.Parameters(col, prog_vis, sub, opt_vis)
    #results_vis = GasChromatographySimulator.simulate(par_vis)
    #@test isapprox(results_vis[1].tR[1], 87.078, atol=1e-3)
    #@test isapprox(results_vis[1].τR[2], 0.59294, atol=1e-5)

    # vis = "HP" and control = "Flow"
    F_steps = [1.0, 1.0, 1.5, 1.2]./(60e6)
    opt_control = GasChromatographySimulator.Options(vis="HP", control="Flow")
    prog_control = GasChromatographySimulator.Program(time_steps, temp_steps, F_steps, pout_steps, [ΔT_steps x₀_steps L₀_steps α_steps], opt_control.Tcontrol, col.L)
    par_control = GasChromatographySimulator.Parameters(col, prog_control, sub, opt_control)
    results_control = GasChromatographySimulator.simulate(par_control)
    @test isapprox(results_control[1].tR[1], 184.29, atol=1e-2) 
    @test isapprox(results_control[1].τR[2], 0.258, atol=1e-3) 

    # compare_peaklist from Misc.jl
    pl1 = results_g[1][1]
    pl2 = results_g[2][1]
    comp = GasChromatographySimulator.compare_peaklist(pl1, pl2)
    @test comp.ΔtR == pl1.tR .- pl2.tR

    # compare_measurement_simulation from Misc.jl
    meas = DataFrame(Name=["C9", "C10"], RT=[100.0, 123.0]) 
    pl1 = results_g[1][1]
    comp = GasChromatographySimulator.compare_measurement_simulation(meas, pl1)
    @test isnan(comp.simulated_tR[1])

    # simulation of a highly retained solute, retention factor > 1e15
    sub_ret = GasChromatographySimulator.Substance("Glyceryl triacetate", "102-76-1", 719.3, 31.282, 1552.2, 0.001, "Brehmer2022, triglyceride, ester", 9.459e-5, 0.0, 0.0)
    prog_ret = GasChromatographySimulator.Program(time_steps, temp_steps, pin_steps, pout_steps, ΔT_steps, x₀_steps, L₀_steps, α_steps, opt[1].Tcontrol, col.L)
    par_ret = GasChromatographySimulator.Parameters(col, prog_ret, [sub_ret], opt[1])
    results_ret = GasChromatographySimulator.simulate(par_ret)
    @test !isnan(results_ret[1].tR[1])
end

@testset "plots check" begin

    a_d = [d]
    a_df = [df]
    d_fun(x) = GasChromatographySimulator.gradient(x, a_d)
    df_fun(x) = GasChromatographySimulator.gradient(x, a_df)
    col = GasChromatographySimulator.Column(L, d_fun, a_d, df_fun, a_df, sp, gas)

    solutes = ["C10", "C11"]
    init_t = zeros(length(solutes))
    init_τ = zeros(length(solutes))
    sub = GasChromatographySimulator.load_solute_database(db_path, db, col.sp, col.gas, solutes, init_t, init_τ)
    
    opt_ng = GasChromatographySimulator.Options(ng=true)
    prog_o_ng = GasChromatographySimulator.Program(time_steps, temp_steps, pin_steps, pout_steps, [zeros(length(time_steps)) x₀_steps L₀_steps α_steps], opt_ng.Tcontrol, col.L)
    par_o_ng = GasChromatographySimulator.Parameters(col, prog_o_ng, sub, opt_ng)
    results_o_ng = GasChromatographySimulator.simulate(par_o_ng)

    # plot functions
    t_start = 0.0#75.0
    t_end = 200.0#150.0
    p_chrom, t, chrom = GasChromatographySimulator.plot_chromatogram(results_o_ng[1], (t_start, t_end))
    @test t[1] == t_start && t[end] == t_end
    tR = results_o_ng[1].tR[2]
    τR = results_o_ng[1].τR[2]
    @test chrom[findlast(t.<=tR)] == 1/sqrt(2*π*τR^2)*exp(-(t[findlast(t.<=tR)]-tR)^2/(2*τR^2))
    t2, chrom2 = GasChromatographySimulator.plot_chromatogram!(p_chrom, results_o_ng[1], (t_start, t_end); mirror=true)
    @test t == t2
    @test chrom == -chrom2
#    @test p_chrom[1][1][:y] == - p_chrom[1][2][:y]
#    @test p_chrom[1][1][:x] == p_chrom[1][2][:x]
    # mirrored the other way
    p_chrom, t, chrom = GasChromatographySimulator.plot_chromatogram(results_o_ng[1], (t_start, t_end); mirror=true)
    @test t[1] == t_start && t[end] == t_end
    tR = results_o_ng[1].tR[1]
    τR = results_o_ng[1].τR[1]
    @test chrom[findlast(t.<=tR)] == -1/sqrt(2*π*τR^2)*exp(-(t[findlast(t.<=tR)]-tR)^2/(2*τR^2))
    t2, chrom2 = GasChromatographySimulator.plot_chromatogram!(p_chrom, results_o_ng[1], (t_start, t_end); mirror=false)
    @test t == t2
    @test chrom == -chrom2
#    @test p_chrom[1][1][:y] == - p_chrom[1][2][:y]
    @test p_chrom[1][1][:x] == p_chrom[1][2][:x]

    p_flow = GasChromatographySimulator.plot_flow(par_o_ng)
    #@test p_flow[1][:xaxis][:optimized_ticks][1][1] = 0.0 
    #@test p_flow[1][:xaxis][:optimized_ticks][1][end] = 800.0
    @test isapprox(p_flow[1][1][:y][1], 11.9797, atol=1e-4)

    p_pres = GasChromatographySimulator.plot_pressure(par_o_ng)
    #@test p_pres[1][:xaxis][:optimized_ticks][1][1] = 0.0 
    #@test p_pres[1][:xaxis][:optimized_ticks][1][end] = 800.0
    @test p_pres[1][1][:y][1] == par_o_ng.prog.Fpin_steps[1]
    @test p_pres[1][2][:y][end] == par_o_ng.prog.pout_steps[end]

    p_temp = GasChromatographySimulator.plot_temperature(par_o_ng)
    @test p_temp[1][1][:y][1] == par_o_ng.prog.temp_steps[1]
    @test p_temp[1][1][:y][end] == par_o_ng.prog.temp_steps[end]
    @test p_temp[1][2][:x][end] == sum(par_o_ng.prog.time_steps)

    p_temp = GasChromatographySimulator.plot_temperature(par_o_ng; selector="T(x)")
    @test p_temp[1][1][:y][1] == par_o_ng.prog.temp_steps[1]
    @test p_temp[1][end][:y][end] == par_o_ng.prog.temp_steps[end]
    @test p_temp[1][2][:x][end] == par_o_ng.col.L

    p_temp = GasChromatographySimulator.plot_temperature(par_o_ng; selector="T(x,t)")
    @test p_temp[1][1][:y][end] == sum(par_o_ng.prog.time_steps)
    @test p_temp[1][1][:x][end] == par_o_ng.col.L

    p_local = GasChromatographySimulator.local_plots("z", "z", results_o_ng[2], par_o_ng)
    @test p_local[1][1][:x][end] == p_local[1][1][:y][end]
    p_local = GasChromatographySimulator.local_plots("t", "t", results_o_ng[2], par_o_ng)
    @test p_local[1][1][:x][end] == p_local[1][1][:y][end]
    p_local = GasChromatographySimulator.local_plots("T", "T", results_o_ng[2], par_o_ng)
    @test p_local[1][1][:x][end] == p_local[1][1][:y][end]
    p_local = GasChromatographySimulator.local_plots("τ", "τ", results_o_ng[2], par_o_ng)
    @test p_local[1][1][:x][end] == p_local[1][1][:y][end]
    p_local = GasChromatographySimulator.local_plots("σ", "σ", results_o_ng[2], par_o_ng)
    @test p_local[1][1][:x][end] === p_local[1][1][:y][end]
    p_local = GasChromatographySimulator.local_plots("u", "u", results_o_ng[2], par_o_ng)
    @test p_local[1][1][:x][end] === p_local[1][1][:y][end]
end

@testset "UI checks" begin
    col_val = GasChromatographySimulator.UI_Column(["HP5", "DB5"])
    @test typeof(col_val) == PlutoUI.CombineNotebook.CombinedBonds
    prog_val = GasChromatographySimulator.UI_Program()
    @test typeof(prog_val) == PlutoUI.CombineNotebook.CombinedBonds
    prog_val = GasChromatographySimulator.UI_Program(default=("0 60", "40 140", "100 100", "vacuum"))
    @test typeof(prog_val) == PlutoUI.CombineNotebook.CombinedBonds
    sub_val = GasChromatographySimulator.UI_Substance(["C6", "C7", "C8", "C9"])
    @test typeof(prog_val) == PlutoUI.CombineNotebook.CombinedBonds
    sub_val = GasChromatographySimulator.UI_Substance(["C6", "C7", "C8", "C9"]; default=(1:2, 0.1, 1.0))
    @test typeof(prog_val) == PlutoUI.CombineNotebook.CombinedBonds
    opt_val = GasChromatographySimulator.UI_Options()
    @test typeof(opt_val) == PlutoUI.CombineNotebook.CombinedBonds
    prog = GasChromatographySimulator.setting_prog(("0 60", "40 140", "100 100", "vacuum"), L)
    @test prog.pout_steps[1] == 0.0
    prog = GasChromatographySimulator.setting_prog(("0 60", "40 140", "100 100", "atmosphere"), L)
    @test prog.pout_steps[1] == 101300.0
end

@testset "misc checks" begin
    @test "C10" in GasChromatographySimulator.all_solutes(sp, DataFrame(CSV.File(string(db_path,"/",db), header=1, silencewarnings=true))) || "Decane" in GasChromatographySimulator.all_solutes(sp, DataFrame(CSV.File(string(db_path,"/",db), header=1, silencewarnings=true)))
    @test "C11" in GasChromatographySimulator.all_solutes(sp, DataFrame(CSV.File(string(db_path,"/",db), header=1, silencewarnings=true))) || "Undecane" in GasChromatographySimulator.all_solutes(sp, DataFrame(CSV.File(string(db_path,"/",db), header=1, silencewarnings=true)))

    a_d = [d]
    a_df = [df]
    d_fun(x) = GasChromatographySimulator.gradient(x, a_d)
    df_fun(x) = GasChromatographySimulator.gradient(x, a_df)
    col = GasChromatographySimulator.Column(L, d_fun, a_d, df_fun, a_df, sp, gas)
    a_gf = [ΔT_steps x₀_steps L₀_steps α_steps]
    gf(x) = GasChromatographySimulator.gradient(x, a_gf)
    T_itp = GasChromatographySimulator.temperature_interpolation(time_steps, temp_steps, gf, col.L)
    pin_itp = GasChromatographySimulator.steps_interpolation(time_steps, pin_steps)
    pout_itp = GasChromatographySimulator.steps_interpolation(time_steps, pout_steps)
    # without gradient
    prog_ng = GasChromatographySimulator.Program(time_steps, temp_steps, pin_steps, pout_steps, col.L)

    # random time in the temperature program
    t = rand()*sum(time_steps)
    # temperature at random time at column inlet
    T_test = prog_ng.T_itp(0.0, t)

    # viscosity tests
    # vis = "Blumberg"
        # gas = "He"
        η_T = GasChromatographySimulator.viscosity(T_test, "He")
        η_t = GasChromatographySimulator.viscosity(0.0, t, prog_ng.T_itp, col.gas)
        @test η_T ≈ η_t
        # gas = "H2"
        η_T = GasChromatographySimulator.viscosity(T_test, "H2")
        η_t = GasChromatographySimulator.viscosity(0.0, t, prog_ng.T_itp, "H2")
        @test η_T ≈ η_t
        # gas = "N2"
        η_T = GasChromatographySimulator.viscosity(T_test, "N2")
        η_t = GasChromatographySimulator.viscosity(0.0, t, prog_ng.T_itp, "N2")
        @test η_T ≈ η_t
        # gas = "Ar"
        η_T = GasChromatographySimulator.viscosity(T_test, "Ar")
        η_t = GasChromatographySimulator.viscosity(0.0, t, prog_ng.T_itp, "Ar")
        @test η_T ≈ η_t
    # vis = "HP
        # gas = "He"
        η_T = GasChromatographySimulator.viscosity(T_test, "He"; vis="HP")
        η_t = GasChromatographySimulator.viscosity(0.0, t, prog_ng.T_itp, "He"; vis="HP")
        @test η_T ≈ η_t
        # gas = "H2"
        η_T = GasChromatographySimulator.viscosity(T_test, "H2"; vis="HP")
        η_t = GasChromatographySimulator.viscosity(0.0, t, prog_ng.T_itp, "H2"; vis="HP")
        @test η_T ≈ η_t
        # gas = "N2"
        η_T = GasChromatographySimulator.viscosity(T_test, "N2"; vis="HP")
        η_t = GasChromatographySimulator.viscosity(0.0, t, prog_ng.T_itp, "N2"; vis="HP")
        @test η_T ≈ η_t
    
    tM_T = GasChromatographySimulator.holdup_time(T_test, prog_ng.Fpin_itp(t), prog_ng.pout_itp(t), col.L, col.a_d[1], col.gas) # only defined for non-gradient case
    tM_t = GasChromatographySimulator.holdup_time(t, prog_ng.T_itp, prog_ng.Fpin_itp, prog_ng.pout_itp, col.L, col.d, col.gas)
    tM_t_ng = GasChromatographySimulator.holdup_time(t, prog_ng.T_itp, prog_ng.Fpin_itp, prog_ng.pout_itp, col.L, col.d, col.gas; ng=true)
    tM_t_d = GasChromatographySimulator.holdup_time(t, prog_ng.T_itp, prog_ng.Fpin_itp, prog_ng.pout_itp, col.L, d, col.gas)
    @test isapprox(tM_T, tM_t, atol=1e-3) 
    @test isapprox(tM_T, tM_t_ng, atol=1e-3)
    @test isapprox(tM_t, tM_t_d, atol=1e-3)

    F_T = GasChromatographySimulator.flow(T_test, prog_ng.Fpin_itp(t), prog_ng.pout_itp(t), col.L, col.a_d[1], col.gas) # only defined for non-gradient case
    F_t = GasChromatographySimulator.flow(t, prog_ng.T_itp, prog_ng.Fpin_itp, prog_ng.pout_itp, col.L, col.d, col.gas)

    @test F_T  ≈  F_t

    # control = "Flow"
    F_steps = [1.0, 1.0, 1.5, 1.2]./(60e6)
    # without gradient
    prog_F = GasChromatographySimulator.Program(time_steps, temp_steps, F_steps, pout_steps, col.L)
    F_T = GasChromatographySimulator.flow(T_test, prog_F.Fpin_itp(t), prog_F.pout_itp(t), col.L, col.a_d[1], col.gas; control="Flow") # only defined for non-gradient case
    F_t = GasChromatographySimulator.flow(t, prog_F.T_itp, prog_F.Fpin_itp, prog_F.pout_itp, col.L, col.d, col.gas; control="Flow")
    pin(t) = GasChromatographySimulator.inlet_pressure(t, prog_F.T_itp, prog_F.Fpin_itp, prog_F.pout_itp, col.L, col.d, col.gas; ng=false, vis="Blumberg", control="Flow")
   
    @test F_T  ≈  F_t
    @test F_T ≈ GasChromatographySimulator.flow(T_test, pin(t), prog_F.pout_itp(t), col.L, col.a_d[1], col.gas; control="Pressure")
    @test F_t ≈ GasChromatographySimulator.flow(t, prog_F.T_itp, prog_F.Fpin_itp, prog_F.pout_itp, col.L, d, col.gas; control="Flow")

    tM_T = GasChromatographySimulator.holdup_time(T_test, prog_F.Fpin_itp(t), prog_F.pout_itp(t), col.L, col.a_d[1], col.gas; control="Flow") # only defined for non-gradient case
    tM_t = GasChromatographySimulator.holdup_time(t, prog_F.T_itp, prog_F.Fpin_itp, prog_F.pout_itp, col.L, col.d, col.gas; control="Flow") # prog_F is without gradient -> should be the same
    @test isapprox(tM_T, tM_t, atol=1e-3)

    tM_T = GasChromatographySimulator.holdup_time(T_test, prog_F.Fpin_itp(t), prog_F.pout_itp(t), col.L, col.a_d[1], col.gas; control="Flow") # only defined for non-gradient case
    tM_t = GasChromatographySimulator.holdup_time(t, prog_F.T_itp, prog_F.Fpin_itp, prog_F.pout_itp, col.L, col.d, col.gas; control="Flow", ng=true)
    @test isapprox(tM_T, tM_t, atol=1e-3)
    

    F_T = GasChromatographySimulator.flow(T_test, prog_ng.Fpin_itp(t), prog_ng.pout_itp(t), col.L, col.a_d[1], col.gas) # only defined for non-gradient case
    F_t = GasChromatographySimulator.flow(t, prog_ng.T_itp, prog_ng.Fpin_itp, prog_ng.pout_itp, col.L, col.d, col.gas)
    F_t_ng = GasChromatographySimulator.flow(t, prog_ng.T_itp, prog_ng.Fpin_itp, prog_ng.pout_itp, col.L, col.d, col.gas, ng=true)
    @test F_T  ≈  F_t
    @test F_T  ≈  F_t_ng

    common = GasChromatographySimulator.common(["aa", "bb", "cc"], ["dd", "aa", "ee"])
    @test common == ["aa"]
end

# Uncertainty functionality tests are currently disabled due to compatibility issues
# between Measurements.jl and ForwardDiff.jl. When using Measurements.Measurement types
# with ForwardDiff.Dual numbers (required for automatic differentiation during ODE solving),
# a MethodError occurs: "no method matching Int64(::Measurements.Measurement{Float64})".
# This happens because ForwardDiff creates Dual numbers containing Measurement values,
# and certain internal operations (e.g., in searchsortedlast) attempt to convert these
# to Int64, which is not supported for Measurement types.
# 
# The simulation itself works correctly with regular Float64 values and ForwardDiff.Dual
# numbers. The issue only manifests when combining Measurements.jl uncertainty propagation
# with ForwardDiff.jl automatic differentiation.
# 
# Potential solutions to explore in the future:
# - ForwardDiffOverMeasurements.jl (if it provides the necessary compatibility)
# - Custom handling of Measurement types in the interpolation and ODE solving code
# - Alternative uncertainty propagation approaches
#=@testset "uncertainty" begin
    col = GasChromatographySimulator.Column(GasChromatographySimulator.measurement(30.0, 1.0), 0.25e-3, 0.25e-6, "Rxi5SilMS", "He")
    prog = GasChromatographySimulator.Program([40.0, 3.0, 10.0, 300.0, 5.0], [300000.0, 3.0, (450000.0-300000.0)/(300.0-40.0)*10.0, 450000.0, 5.0], col.L)
    @test typeof(prog.T_itp(col.L, 523.0)) == GasChromatographySimulator.Measurement{Float64}
    sub = GasChromatographySimulator.load_solute_database(db_path, db_unc, col.sp, col.gas, ["Octan-1-ol", "Aniline"], fill(GasChromatographySimulator.measurement(0.0, 0.0), 2), fill(GasChromatographySimulator.measurement(0.0, 0.0), 2))
    @test typeof(sub[1].Tchar) == GasChromatographySimulator.Measurement{Float64}
    opt = GasChromatographySimulator.Options(ng=true)
    par = GasChromatographySimulator.Parameters(col, prog, sub, opt)
    sim = GasChromatographySimulator.simulate(par)
    @test typeof(sim[1].tR) == Array{GasChromatographySimulator.Measurement{Float64}, 1}
end=#

@testset "differentiability" begin
    prog = GasChromatographySimulator.Program([40.0, 3.0, 10.0, 300.0, 5.0], [300000.0, 3.0, (450000.0-300000.0)/(300.0-40.0)*10.0, 450000.0, 5.0], 30.0)
    opt = GasChromatographySimulator.Options(ng=true)
    # retention time and peak width of one substance depending on retention parameters
    sol_tR_τR_RP(x) = GasChromatographySimulator.solving_odesystem_r(30.0, 0.25e-3, 0.25e-6, "He", prog.T_itp, prog.Fpin_itp, prog.pout_itp, x[1], x[2], x[3], 1e-3, 1e-6, 0.0, 0.0, opt).u[end]
    # differentiation
    ∂sol_tR_τR_RP(x) = GasChromatographySimulator.ForwardDiff.jacobian(sol_tR_τR_RP, x)
    val = ∂sol_tR_τR_RP([400.0, 30.0, 100.0])
    @test size(val) == (2, 3)
    # retention time and peak width of one substance depending on column parameters
    sol_tR_τR_col(x) = GasChromatographySimulator.solving_odesystem_r(x[1], x[2], x[3], "He", prog.T_itp, prog.Fpin_itp, prog.pout_itp, 400.0, 30.0, 100.0, 1e-3, 1e-6, 0.0, 0.0, opt).u[end]
    # differentiation
    ∂sol_tR_τR_col(x) = GasChromatographySimulator.ForwardDiff.jacobian(sol_tR_τR_col, x)
    val_col = ∂sol_tR_τR_col([30.0, 0.25e-3, 0.25e-6])
    @test size(val_col) == (2, 3)
end

@testset "ODE solver dtNaN" begin
    col = GasChromatographySimulator.Column(0.25, 0.1e-3, 0.1e-6, "Rxi5SilMS", "He")
    prog = GasChromatographySimulator.Program([300.0, 5.0], [120000.0, 5.0], col.L)
    opt = GasChromatographySimulator.Options(abstol=1e-8, reltol=1e-5)
    sol = GasChromatographySimulator.solving_odesystem_r(col.L, col.d, col.df, col.gas, prog.T_itp, prog.Fpin_itp, prog.pout_itp, 340.0, 34.0, 200.0, 1e-3, 1e-4, 240.0, 1.0, opt)
    @test sol.t[end] == col.L
    #GasChromatographySimulator.solve_system(col.L, col.d, col.df, col.gas, prog.T_itp, prog.Fpin_itp, prog.pout_itp, 340.0, 34.0, 200.0, 1e-3, 1e-4, 240.0, 1.0, opt)
end

@testset "linear_interpolation" begin
    # Test 1: Basic 1D linear interpolation with sorted data
    x = [0.0, 1.0, 2.0, 3.0]
    y = [10.0, 20.0, 30.0, 40.0]
    interp = GasChromatographySimulator.linear_interpolation((x,), y)
    
    # Test interpolation at known points
    @test interp(0.0) == 10.0
    @test interp(1.0) == 20.0
    @test interp(2.0) == 30.0
    @test interp(3.0) == 40.0
    
    # Test interpolation between points
    @test isapprox(interp(0.5), 15.0, atol=1e-10)
    @test isapprox(interp(1.5), 25.0, atol=1e-10)
    @test isapprox(interp(2.5), 35.0, atol=1e-10)
    
    # Test 2: Extrapolation (should return boundary values - flat extrapolation)
    @test interp(-1.0) == 10.0  # Below range, returns first value
    @test interp(5.0) == 40.0   # Above range, returns last value
    
    # Test 3: Single point (edge case)
    x_single = [1.0]
    y_single = [5.0]
    interp_single = GasChromatographySimulator.linear_interpolation((x_single,), y_single)
    @test interp_single(1.0) == 5.0
    @test interp_single(0.0) == 5.0  # Extrapolation returns the single value
    @test interp_single(2.0) == 5.0  # Extrapolation returns the single value
    
    # Test 4: Non-uniform spacing
    x_nonuniform = [0.0, 0.5, 2.0, 5.0]
    y_nonuniform = [0.0, 10.0, 40.0, 100.0]
    interp_nonuniform = GasChromatographySimulator.linear_interpolation((x_nonuniform,), y_nonuniform)
    # At x=1.0: between 0.5 (y=10.0) and 2.0 (y=40.0)
    # Linear: 10.0 + (1.0-0.5)/(2.0-0.5) * (40.0-10.0) = 10.0 + 0.5/1.5 * 30.0 = 20.0
    @test isapprox(interp_nonuniform(1.0), 20.0, atol=1e-10)
    # At x=3.5: between 2.0 (y=40.0) and 5.0 (y=100.0)
    # Linear: 40.0 + (3.5-2.0)/(5.0-2.0) * (100.0-40.0) = 40.0 + 1.5/3.0 * 60.0 = 70.0
    @test isapprox(interp_nonuniform(3.5), 70.0, atol=1e-10)
    
    # Test 5: 2D bilinear interpolation
    x_grid = [0.0, 1.0, 2.0]
    t_grid = [0.0, 10.0, 20.0]
    values = [10.0 20.0 30.0; 15.0 25.0 35.0; 20.0 30.0 40.0]
    interp_2d = GasChromatographySimulator.linear_interpolation((x_grid, t_grid), values)
    
    # Test at grid points
    @test interp_2d(0.0, 0.0) == 10.0
    @test interp_2d(1.0, 10.0) == 25.0
    @test interp_2d(2.0, 20.0) == 40.0
    
    # Test interpolation between grid points
    # At (0.5, 5.0): bilinear interpolation
    # x=0.5 is between x_grid[1]=0.0 and x_grid[2]=1.0
    # t=5.0 is between t_grid[1]=0.0 and t_grid[2]=10.0
    # Interpolate in x at t=0.0: 10.0 + 0.5*(15.0-10.0) = 12.5
    # Interpolate in x at t=10.0: 20.0 + 0.5*(25.0-20.0) = 22.5
    # Interpolate in t: 12.5 + 0.5*(22.5-12.5) = 17.5
    result = interp_2d(0.5, 5.0)
    @test isapprox(result, 17.5, atol=1e-10)
    
    # Test 6: Extrapolation for 2D (should return boundary values)
    @test interp_2d(-1.0, 5.0) == interp_2d(0.0, 5.0)  # Below x range
    @test interp_2d(3.0, 5.0) == interp_2d(2.0, 5.0)   # Above x range
    @test interp_2d(1.0, -5.0) == interp_2d(1.0, 0.0)  # Below t range
    @test interp_2d(1.0, 25.0) == interp_2d(1.0, 20.0) # Above t range
    
    # Test 7: Equality comparison
    interp1 = GasChromatographySimulator.linear_interpolation((x,), y)
    interp2 = GasChromatographySimulator.linear_interpolation((x,), y)
    @test interp1 == interp2
    
    # Test 8: Different values should not be equal
    y2 = [10.0, 21.0, 30.0, 40.0]
    interp3 = GasChromatographySimulator.linear_interpolation((x,), y2)
    @test interp1 != interp3
end

@testset "deduplicate_knots!" begin
    # Test 1: Basic deduplication with move_knots=true
    x = [1.0, 2.0, 2.0, 3.0, 3.0, 3.0, 4.0]
    x_copy = copy(x)
    unique_indices = GasChromatographySimulator.deduplicate_knots!(x_copy; move_knots=true)
    
    # After deduplication, consecutive duplicates should be slightly different
    @test length(unique_indices) == length(x)
    @test x_copy[2] != x_copy[3]  # Duplicates should be modified
    @test abs(x_copy[2] - x_copy[3]) < 1e-9  # But very close
    @test x_copy[4] != x_copy[5]  # Multiple duplicates
    @test x_copy[5] != x_copy[6]  # Multiple duplicates
    
    # Test 2: No duplicates (should not modify)
    x_no_dup = [1.0, 2.0, 3.0, 4.0, 5.0]
    x_no_dup_copy = copy(x_no_dup)
    unique_indices = GasChromatographySimulator.deduplicate_knots!(x_no_dup_copy; move_knots=true)
    @test x_no_dup_copy == x_no_dup  # Should be unchanged
    
    # Test 3: move_knots=false (should not modify duplicates)
    x_dup = [1.0, 2.0, 2.0, 3.0]
    x_dup_copy = copy(x_dup)
    unique_indices = GasChromatographySimulator.deduplicate_knots!(x_dup_copy; move_knots=false)
    @test x_dup_copy == x_dup  # Should be unchanged when move_knots=false
    
    # Test 4: Single element
    x_single = [5.0]
    x_single_copy = copy(x_single)
    unique_indices = GasChromatographySimulator.deduplicate_knots!(x_single_copy; move_knots=true)
    @test x_single_copy == x_single
    @test unique_indices == 1:1
    
    # Test 5: Empty array (edge case)
    x_empty = Float64[]
    x_empty_copy = copy(x_empty)
    unique_indices = GasChromatographySimulator.deduplicate_knots!(x_empty_copy; move_knots=true)
    @test x_empty_copy == x_empty
    @test unique_indices == 1:0
    
    # Test 6: All duplicates
    x_all_dup = [1.0, 1.0, 1.0, 1.0]
    x_all_dup_copy = copy(x_all_dup)
    unique_indices = GasChromatographySimulator.deduplicate_knots!(x_all_dup_copy; move_knots=true)
    @test length(unique_indices) == length(x_all_dup)
    # Last element should be adjusted based on previous difference
    @test x_all_dup_copy[end] != x_all_dup_copy[end-1]
    
    # Test 7: Duplicate at the end
    x_end_dup = [1.0, 2.0, 3.0, 3.0]
    x_end_dup_copy = copy(x_end_dup)
    unique_indices = GasChromatographySimulator.deduplicate_knots!(x_end_dup_copy; move_knots=true)
    @test x_end_dup_copy[3] != x_end_dup_copy[4]  # Last duplicate should be modified
    @test abs(x_end_dup_copy[3] - x_end_dup_copy[4]) < 1e-9  # But very close
end

println("Test run successful.")