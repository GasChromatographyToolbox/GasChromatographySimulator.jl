#=using GasChromatographySimulator
using CSV, DataFrames
using Plots

db_path = "/Users/janleppert/Documents/GitHub/GasChromatographySimulator/data"#"/../../data"
db_name = "Database_Leppert2020b.csv"

# conventional GC
opt = GasChromatographySimulator.Options(abstol=1e-8, reltol=1e-5, ng=false) # setting ng=true for conventional GC, use of simplified model, slightly faster than ng=false
col = GasChromatographySimulator.Column(11.18, 0.104e-3, 0.104e-6, "FS5ms", "H2")

# solutes
db_dataframe = DataFrame(CSV.File(string(db_path,"/",db_name), header=1, silencewarnings=true)) # load the database
solutes = GasChromatographySimulator.all_solutes(col.sp, db_dataframe) # all solutes of the database for this stationary phase 
sub = GasChromatographySimulator.load_solute_database(db_dataframe, col.sp, col.gas, solutes, zeros(length(solutes)), zeros(length(solutes))) # load all solutes, injection at start of the program, ideal injection (what was used in the paper?)

# programs
prog = Array{GasChromatographySimulator.Program}(undef, 4) # 4 different programs
prog[1] = GasChromatographySimulator.Program([0.0, 60.0, 1560.0, 180.0], [40.0, 40.0, 300.0, 300.0], 411564.0*ones(4), 101300.0.*ones(4), col.L)
prog[2] = GasChromatographySimulator.Program([0.0, 60.0, 1040.0, 180.0], [40.0, 40.0, 300.0, 300.0], 411564.0*ones(4), 101300.0.*ones(4), col.L)
prog[3] = GasChromatographySimulator.Program([0.0, 60.0, 3120.0, 180.0], [40.0, 40.0, 300.0, 300.0], 411564.0*ones(4), 101300.0.*ones(4), col.L)
prog[4] = GasChromatographySimulator.Program([0.0, 60.0, 1680.0, 60.0, 360.0, 60.0], [40.0, 40.0, 180.0, 180.0, 300.0, 300.0], 411564.0*ones(6), 101300.0.*ones(6), col.L)

# parameters
par = Array{GasChromatographySimulator.Parameters}(undef, 4)
for i=1:4
    par[i] = GasChromatographySimulator.Parameters(col, prog[i], sub, opt)
end

peaklist = Array{DataFrame}(undef, 4)
sol = Array{Any}(undef, 4)

for i=1:4
    peaklist[i], sol[i] = GasChromatographySimulator.simulate(par[i])
end

# load measured data:
data_path = "/Users/janleppert/Documents/GitHub/GasChromatographySimulator/data/measurements"#"/../../data"
data_name = ["Leppert2020b_measured_RT_progA.csv", "Leppert2020b_measured_RT_progB.csv", "Leppert2020b_measured_RT_progC.csv", "Leppert2020b_measured_RT_progD.csv"]

meas_data = Array{DataFrame}(undef, 4)
for i=1:4
    meas_data[i] = DataFrame(CSV.File(string(data_path,"/",data_name[i]), header=1, silencewarnings=true)) # load the measured retention times
    meas_data[i][!, 2] = meas_data[i][!, 2] .* 60.0 # conversion from min -> s
    rename!(meas_data[i], [:Name, :tR, :τR])
end

compare_ = Array{DataFrame}(undef, 4)
for i=1:4
    compare_[i] = GasChromatographySimulator.compare_peaklist(meas_data[i], peaklist[i])
end

data_chrom = "Leppert2020b_measured_Chrom_progD.csv"

meas_chrom = DataFrame(CSV.File(string(data_path, "/", data_chrom), header=1, silencewarnings=true))

p_chrom, t, chrom = GasChromatographySimulator.plot_chromatogram(peaklist[4], (0.0, 2200.0); annotation=false, number=true, mirror=true, offset=0.0)

p_chrom
using Plots
plot!(p_chrom,meas_chrom[!,1], meas_chrom[!,2].*400.0.+0.1)
ylims!(-1.6,1.6)
xlims!(0.0,round(meas_chrom[end,1];sigdigits=2))


# with gradient
opt = GasChromatographySimulator.Options() # setting ng=true for conventional GC, use of simplified model, slightly faster than ng=false
col = GasChromatographySimulator.Column(2.05, 0.104e-3, 0.104e-6, "FS5ms", "He")
    


prog_settings = DataFrame(CSV.File(string(data_path,"/x90.csv"), header=1))
time_step = prog_settings.Deltat
temp_step = prog_settings.T
ΔT_steps = prog_settings.DeltaT
pin_steps = prog_settings.pinj.*1000.0 .+ 101300.0
pout_steps = prog_settings.pdet.*1000.0
α_steps = -3.0.*ones(length(ΔT_steps))
x₀_steps = zeros(length(ΔT_steps))
L₀_steps = col.L.*ones(length(ΔT_steps))
prog_high_grad = GasChromatographySimulator.Program(time_step, temp_step, pin_steps, pout_steps, ΔT_steps, x₀_steps, L₀_steps, α_steps, "outlet", col.L)


db_dataframe = DataFrame(CSV.File(string(db_path,"/",db_name), header=1, silencewarnings=true)) # load the database
solutes = GasChromatographySimulator.all_solutes(col.sp, db_dataframe) # all solutes of the database for this stationary phase 
sub = GasChromatographySimulator.load_solute_database(db_dataframe, col.sp, col.gas, solutes, zeros(length(solutes)), zeros(length(solutes))) # load all solutes, injection at start of the program, ideal injection (what was used in the paper?)

par = GasChromatographySimulator.Parameters(col, prog_high_grad, sub, opt)

p_flow_tg = GasChromatographySimulator.plot_flow(par)
p_press_tg = GasChromatographySimulator.plot_pressure(par)
p_temp_tg = GasChromatographySimulator.plot_temperature(par)
plot!(p_temp_tg, label=:left)
l = @layout([a{0.65w} [b; c]])
p_TpF_tg = plot(p_temp_tg, p_press_tg, p_flow_tg, layout=l)


peaklist, sol = GasChromatographySimulator.simulate(par);

peaklist


meas_data = DataFrame(CSV.File(string(data_path,"/Leppert2020b_measured_RT_med_gradient.csv"), header=1, silencewarnings=true))
meas_data[!, 3] = meas_data[!, 3] ./ 1000.0 # conversion from ms -> s
rename!(meas_data, [:Name, :tR, :τR])
compare = GasChromatographySimulator.compare_peaklist(meas_data, peaklist)

data_chrom = "Leppert2020b_measured_Chrom_med_gradient_x90.csv"

meas_chrom = DataFrame(CSV.File(string(data_path, "/", data_chrom), header=1, silencewarnings=true))


p_chrom, t, chrom = GasChromatographySimulator.plot_chromatogram(peaklist, (0.0, 70.0); annotation=false, number=true, mirror=true, offset=0.0)

plot!(p_chrom,meas_chrom[!,1].*60.0, meas_chrom[!,2].*8e-5)
ylims!(-13,13)
xlims!(0.0,55.0)
=#




# ------------
## Simulation of measurements

using GasChromatographySimulator # hide
using DataFrames, CSV # hide
using Plots # hide
opt = GasChromatographySimulator.Options(ng=true)

col = GasChromatographySimulator.Column(11.18, 0.104e-3, 0.104e-6, "FS5ms", "H2")

prog_D = GasChromatographySimulator.Program([0.0, 60.0, 1680.0, 60.0, 360.0, 60.0], [40.0, 40.0, 180.0, 180.0, 300.0, 300.0], 411564.0*ones(6), 101300.0.*ones(6), col.L)

db_dataframe = DataFrame(CSV.File("/Users/janleppert/Documents/GitHub/GasChromatographySimulator/data/Leppert2020b/Database_Leppert2020b.csv", header=1, silencewarnings=true))

solutes = GasChromatographySimulator.all_solutes(col.sp, db_dataframe)

t₀ = zeros(length(solutes))
τ₀ = zeros(length(solutes))
sub = GasChromatographySimulator.load_solute_database(db_dataframe, col.sp, col.gas, solutes, t₀, τ₀)

par = GasChromatographySimulator.Parameters(col, prog_D, sub, opt)

p_flow = GasChromatographySimulator.plot_flow(par)
p_press = GasChromatographySimulator.plot_pressure(par)
p_temp = GasChromatographySimulator.plot_temperature(par)
l = @layout([a{0.65w} [b; c]])
p_TpF = plot(p_temp, p_press, p_flow, layout=l)

peaklist, sol = GasChromatographySimulator.simulate(par);
peaklist

measurement_D = DataFrame(CSV.File("/Users/janleppert/Documents/GitHub/GasChromatographySimulator/data/Leppert2020b/Leppert2020b_measured_RT_progD.csv", header=1, silencewarnings=true))
measurement_D[!, 2] = measurement_D[!, 2] .* 60.0 # conversion from min -> s
rename!(measurement_D, [:Name, :tR, :τR])
compare = GasChromatographySimulator.compare_peaklist(measurement_D, peaklist)

gr()
chrom_D = DataFrame(CSV.File("/Users/janleppert/Documents/GitHub/GasChromatographySimulator/data/Leppert2020b/Leppert2020b_measured_Chrom_progD.csv", header=1, silencewarnings=true))
p_chrom, t, chrom = GasChromatographySimulator.plot_chromatogram(peaklist, (0.0, round(chrom_D[end,1];sigdigits=2)); annotation=false, number=true, mirror=true, offset=0.0)
plot!(p_chrom, chrom_D[!,1], chrom_D[!,2].*400.0.+0.1)
ylims!(-1.6,1.6)
xlims!(0.0,round(chrom_D[end,1];sigdigits=2))
p_chrom
#savefig(p_chrom,"compare_conventionalGC.svg")

opt_tg = GasChromatographySimulator.Options()

col_tg = GasChromatographySimulator.Column(2.05, 0.104e-3, 0.104e-6, "FS5ms", "He")

prog_settings = DataFrame(CSV.File("/Users/janleppert/Documents/GitHub/GasChromatographySimulator/data/Leppert2020b/Leppert2020b_prog_settings_med_gradient_x90.csv", header=1, silencewarnings=true))
# reduce settings to every 20th point 
time = cumsum(prog_settings.Deltat)[1:20:end]
time_steps = Array{Float64}(undef, length(time))
for i=2:length(time)
    time_steps[i] = time[i]-time[i-1]
end
time_steps[1] = 0.0

#time_steps = prog_settings.Deltat
temp_steps = prog_settings.T[1:20:end]
ΔT_steps = prog_settings.DeltaT[1:20:end]
pin_steps = prog_settings.pinj[1:20:end].*1000.0 .+ 101300.0
pout_steps = prog_settings.pdet[1:20:end].*1000.0
α_steps = -3.0.*ones(length(ΔT_steps))
x₀_steps = zeros(length(ΔT_steps))
L₀_steps = col_tg.L.*ones(length(ΔT_steps))
prog_med_grad = GasChromatographySimulator.Program(time_steps, temp_steps, pin_steps, pout_steps, ΔT_steps, x₀_steps, L₀_steps, α_steps, "outlet", col_tg.L)

par_tg = GasChromatographySimulator.Parameters(col_tg, prog_med_grad, sub, opt_tg)

plotly()
p_flow_tg = GasChromatographySimulator.plot_flow(par_tg)
p_press_tg = GasChromatographySimulator.plot_pressure(par_tg)
p_temp_tg = GasChromatographySimulator.plot_temperature(par_tg)
l = @layout([a{0.65w} [b; c]])
p_TpF_tg = plot(p_temp_tg, p_press_tg, p_flow_tg, layout=l)



peaklist_tg, sol_tg = GasChromatographySimulator.simulate(par_tg);
peaklist_tg

measurement_tg = DataFrame(CSV.File("/Users/janleppert/Documents/GitHub/GasChromatographySimulator/data/Leppert2020b/Leppert2020b_measured_RT_med_gradient.csv", header=1, silencewarnings=true))
measurement_tg[!, 3] = measurement_tg[!, 3] ./ 1000.0 # conversion from ms -> s
rename!(measurement_tg, [:Name, :tR, :τR])

compare_tg = GasChromatographySimulator.compare_peaklist(measurement_tg, peaklist_tg)

gr()
chrom_tg = DataFrame(CSV.File("/Users/janleppert/Documents/GitHub/GasChromatographySimulator/data/Leppert2020b/Leppert2020b_measured_Chrom_med_gradient_x90.csv", header=1, silencewarnings=true))
p_chrom_tg, t_, chrom_ = GasChromatographySimulator.plot_chromatogram(peaklist_tg, (0.0, 55.0); annotation=false, number=true, mirror=true, offset=0.0)
plot!(p_chrom_tg, chrom_tg[!,1].*60.0, chrom_tg[!,2].*8e-5)
ylims!(-13,13)
xlims!(0.0,55.0)
p_chrom_tg
#savefig(p_chrom_tg,"compare_TGGC.svg")