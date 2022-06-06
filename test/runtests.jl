using Test, GasChromatographySimulator, PlutoUI

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
    col_c = GasChromatographySimulator.constructor_System(L, a_d[1], a_df[1], sp, gas)
    @test col.d(col.L) == col_c.d(col_c.L)

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
    db_path = string(@__DIR__, "/data")
    db = "Database_test.csv"
    solutes = ["C10", "C11"]
    init_t = zeros(length(solutes))
    init_τ = zeros(length(solutes))
    sub = GasChromatographySimulator.load_solute_database(db_path, db, col.sp, col.gas, solutes, init_t, init_τ)
    @test sub[1].CAS == "124-18-5"
    @test sub[2].Tchar == 124.746 + 273.15

    # no stationary phase & other gas for diffusivity
    # H2
    sub_0 = GasChromatographySimulator.load_solute_database(db_path, db, "", "H2", solutes, init_t, init_τ)
    @test sub_0[1].Tchar == 1.0
    @test sub_0[2].ann == "no sp"
    @test round(sub_0[1].Dag; sigdigits=5) == 0.00011885
    # N2
    sub_0 = GasChromatographySimulator.load_solute_database(db_path, db, "", "N2", solutes, init_t, init_τ)
    @test round(sub_0[1].Dag; sigdigits=5) == 2.8398e-5
    # Ar
    sub_0 = GasChromatographySimulator.load_solute_database(db_path, db, "", "Ar", solutes, init_t, init_τ)
    @test round(sub_0[1].Dag; sigdigits=5) == 2.5268e-5
    # test for error-cases of GasChromatographySimulator.load_solute_database

    # test for new database format
    db_new = "Database_test_new_format.csv"
    sub_new = GasChromatographySimulator.load_solute_database(db_path, db_new, "Wax", "He", ["C14", "Decyl acetate", "Hexadecane", "C15", "Methyl myristate"], [1.0, 2.0, 3.0, 4.0, 5.0], [0.1, 0.2, 0.3, 0.4, 0.5])
    sub_old = GasChromatographySimulator.load_solute_database(db_path, db, "Wax", "He", ["C14", "Decyl acetate", "Hexadecane", "C15", "Methyl myristate"], [1.0, 2.0, 3.0, 4.0, 5.0], [0.1, 0.2, 0.3, 0.4, 0.5])
    @test sub_new[1].CAS == sub_old[1].CAS
    @test sub_new[2].t₀ == 4.0

    # parameters
    par = GasChromatographySimulator.Parameters(col, prog, sub, opt)
    @test par.col.L == L
    @test par.prog.gf(0.0) == par.prog.gf(L) + par.prog.a_gf[:,1]
    @test par.prog.T_itp(0.0, 0.0) == temp_steps[1] + 273.15
    @test par.prog.T_itp(L, sum(time_steps)) == temp_steps[end] - ΔT_steps[end] + 273.15
    @test par.sub[1].Dag == GasChromatographySimulator.diffusivity(par.sub[1].CAS, "He")
    @test  isapprox(par.sub[2].Dag, GasChromatographySimulator.diffusivity(156.31, 11, 24, 0, 0, 0, "He"), atol=1e-6)
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
    a_d = [d]
    a_df = [df]
    d_fun(x) = GasChromatographySimulator.gradient(x, a_d)
    df_fun(x) = GasChromatographySimulator.gradient(x, a_df)
    col = GasChromatographySimulator.Column(L, d_fun, a_d, df_fun, a_df, sp, gas)

    
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
    @test isapprox(results_g[1][1].tR[1], 123.186, atol=1e-3)
    @test isapprox(results_g[2][1].tR[1], 51.4564, atol=1e-4)

    #@test isapprox(results_g[1][1].τR[2], 0.548163, atol=1e-5)
    @test isapprox(results_g[1][1].τR[2], 0.547683, atol=1e-5)
    @test isapprox(results_g[2][1].τR[2], 0.524873, atol=1e-4)

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
    @test isapprox(results_o_ng[1].tR[1], 87.401, atol=1e-3)
    @test isapprox(results_o_ng[1].τR[2], 0.59028, atol=1e-5)

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
    @test isapprox(results_control[1].tR[1], 184.280, atol=1e-3) 
    @test isapprox(results_control[1].τR[2], 0.25771, atol=1e-5) 

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
    @test p_chrom[1][1][:y] == - p_chrom[1][2][:y]
    @test p_chrom[1][1][:x] == p_chrom[1][2][:x]
    # mirrored the other way
    p_chrom, t, chrom = GasChromatographySimulator.plot_chromatogram(results_o_ng[1], (t_start, t_end); mirror=true)
    @test t[1] == t_start && t[end] == t_end
    tR = results_o_ng[1].tR[1]
    τR = results_o_ng[1].τR[1]
    t0 = round(tR; digits=1)
    @test chrom[findfirst(t.==t0)] == -1/sqrt(2*π*τR^2)*exp(-(t0-tR)^2/(2*τR^2))
    t2, chrom2 = GasChromatographySimulator.plot_chromatogram!(p_chrom, results_o_ng[1], (t_start, t_end); mirror=false)
    @test t == t2
    @test chrom == -chrom2
    @test p_chrom[1][1][:y] == - p_chrom[1][2][:y]
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
    @test "C10" in GasChromatographySimulator.all_solutes(sp, DataFrame(CSV.File(string(db_path,"/",db), header=1, silencewarnings=true)))
    @test "C11" in GasChromatographySimulator.all_solutes(sp, DataFrame(CSV.File(string(db_path,"/",db), header=1, silencewarnings=true)))

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
    @test tM_T ≈ tM_t
    @test tM_T ≈ tM_t_ng

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

    tM_T = GasChromatographySimulator.holdup_time(T_test, prog_F.Fpin_itp(t), prog_F.pout_itp(t), col.L, col.a_d[1], col.gas; control="Flow") # only defined for non-gradient case
    tM_t = GasChromatographySimulator.holdup_time(t, prog_F.T_itp, prog_F.Fpin_itp, prog_F.pout_itp, col.L, col.d, col.gas; control="Flow") # prog_F is without gradient -> should be the same
    @test tM_T ≈ tM_t

    tM_T = GasChromatographySimulator.holdup_time(T_test, prog_F.Fpin_itp(t), prog_F.pout_itp(t), col.L, col.a_d[1], col.gas; control="Flow") # only defined for non-gradient case
    tM_t = GasChromatographySimulator.holdup_time(t, prog_F.T_itp, prog_F.Fpin_itp, prog_F.pout_itp, col.L, col.d, col.gas; control="Flow", ng=true)
    @test tM_T ≈ tM_t
    

    F_T = GasChromatographySimulator.flow(T_test, prog_ng.Fpin_itp(t), prog_ng.pout_itp(t), col.L, col.a_d[1], col.gas) # only defined for non-gradient case
    F_t = GasChromatographySimulator.flow(t, prog_ng.T_itp, prog_ng.Fpin_itp, prog_ng.pout_itp, col.L, col.d, col.gas)
    F_t_ng = GasChromatographySimulator.flow(t, prog_ng.T_itp, prog_ng.Fpin_itp, prog_ng.pout_itp, col.L, col.d, col.gas, ng=true)
    @test F_T  ≈  F_t
    @test F_T  ≈  F_t_ng

    common = GasChromatographySimulator.common(["aa", "bb", "cc"], ["dd", "aa", "ee"])
    @test common == ["aa"]
end

println("Test run successful.")