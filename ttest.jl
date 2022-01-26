using Pkg
Pkg.activate("/Users/janleppert/Documents/GitHub/GasChromatographySimulator")
using GasChromatographySimulator
using Plots
# make some graphs for the documentation

sys = GasChromatographySimulator.constructor_System(10.0, 0.25e-3, 0.25e-6, "SLB5ms", "He")
# without gradient
prog_ng = GasChromatographySimulator.Program(   [0.0, 60.0, 300.0, 120.0],
                                                [40.0, 40.0, 300.0, 300.0],
                                                [18.0, 18.0, 98.0, 98.0].*1000.0 .+ 101300.0,
                                                [101.3, 101.3, 101.3, 101.3].*1000.0,
                                                sys.L
                                            )

sub = GasChromatographySimulator.load_solute_database("./data", "Database_test.csv", sys.sp, sys.gas, ["C10", "C11"], [0.0, 0.0], [0.0, 0.0])

opt = GasChromatographySimulator.Options()

par_ng = GasChromatographySimulator.Parameters(sys, prog_ng, sub, opt)

# temperature plot T(t)
t = 0.0:sum(par_ng.prog.time_steps)/1000:sum(par_ng.prog.time_steps)
T0 = par_ng.prog.T_itp.(0.0, t).-273.15
Tplot = plot(t, T0, ylabel="T in °C", xticks=0:60:480, legend=false, label="inlet", c=:red, size=(400,300))
#xgrid!(Tplot, :on, :grey, 2, :solid, 0.8)
#ygrid!(Tplot, :on, :grey, 2, :solid, 0.8)
#TL = par_ng.prog.T_itp.(par_ng.sys.L, nt).-273.15
#plot!(Tplot, nt, TL, label="outlet", c=:blue)

pin = Array{Float64}(undef, length(t))
pout = Array{Float64}(undef, length(t))
for i=1:length(t)
    pin[i] = par_ng.prog.pin_itp(t[i])
    pout[i] = par_ng.prog.pout_itp(t[i])
end
ppress = plot(t, pin./1000.0, xlabel="time in s", xticks=0:60:480, ylabel="pressure in kPa", label="inlet", size=(400,300), legend=false, c=:green)
plot!(ppress, t, pout./1000.0, label="outlet", c=:lime)
#xgrid!(ppress, :on, :grey, 2, :solid, 0.8)
#ygrid!(ppress, :on, :grey, 2, :solid, 0.8)

l = @layout([a; b])
p = plot(Tplot, ppress, layout=l, size=(500, 300))

savefig(p, "prog_nograd.svg")



# with gradient
prog_g = GasChromatographySimulator.Program([0.0, 60.0, 150.0, 150.0, 120.0],
                                            [40.0, 40.0, 170.0, 300.0, 300.0],
                                            150000.0.*ones(5),
                                            zeros(5),
                                            [[0.0, 0.0, 60.0, 60.0, 20.0] zeros(5) sys.L.*ones(5) [0.0, 0.0, -2.0, -5.0, -5.0]],
                                            "inlet",
                                            sys.L)

par_g = GasChromatographySimulator.Parameters(sys, prog_g, sub, opt)

t = 0.0:sum(par_g.prog.time_steps)/1000:sum(par_g.prog.time_steps)
T0 = par_g.prog.T_itp.(0.0, t).-273.15
Tplot_g = plot(t, T0, ylabel="T in °C", xlabel="time in s", xticks=0:60:480, legend=false, label="inlet", c=:red, size=(500,300))
TL = par_g.prog.T_itp.(par_g.sys.L, t).-273.15
plot!(Tplot_g, t, TL, c=:blue)

savefig(Tplot_g, "prog_grad.svg")

# with gradient, outlet
prog_g = GasChromatographySimulator.Program([0.0, 60.0, 150.0, 150.0, 120.0],
                                            [40.0, 40.0, 170.0, 300.0, 300.0],
                                            150000.0.*ones(5),
                                            zeros(5),
                                            [[0.0, 0.0, 60.0, 60.0, 20.0] zeros(5) sys.L.*ones(5) [0.0, 0.0, -2.0, -5.0, -5.0]],
                                            "outlet",
                                            sys.L)

par_g = GasChromatographySimulator.Parameters(sys, prog_g, sub, GasChromatographySimulator.Options(Tcontrol="outlet"))

t = 0.0:sum(par_g.prog.time_steps)/1000:sum(par_g.prog.time_steps)
T0 = par_g.prog.T_itp.(0.0, t).-273.15
Tplot_g = plot(t, T0, ylabel="T in °C", xlabel="time in s", xticks=0:60:480, legend=false, label="inlet", c=:red, size=(500,300))
TL = par_g.prog.T_itp.(par_g.sys.L, t).-273.15
plot!(Tplot_g, t, TL, c=:blue)

savefig(Tplot_g, "prog_grad_outlet.svg")

# different gradient profiles - α
prog_g = GasChromatographySimulator.Program([0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0],
                                            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                            150000.0.*ones(7),
                                            zeros(7),
                                            [50.0.*ones(7) zeros(7) sys.L.*ones(7) [-10.0, -5.0, -2.5, 0.0, 2.5, 5.0, 10.0]],
                                            "inlet",
                                            sys.L)
par_g = GasChromatographySimulator.Parameters(sys, prog_g, sub, GasChromatographySimulator.Options())

x = 0.0:par_g.sys.L/1000:par_g.sys.L
grad = Array{Array{Float64,1}}(undef, 7)
gradplot = plot(xlabel="position in m", ylabel="ΔT in °C")
for i=1:7
    grad[i] = par_g.prog.T_itp.(x, cumsum(par_g.prog.time_steps)[i]).-273.15
    plot!(gradplot, x, grad[i], label="α=$(par_g.prog.a_gf[i,4])")
end
gradplot

savefig(gradplot, "grad_alpha.svg")

# different gradient profiles - x₀
prog_g = GasChromatographySimulator.Program([0.0, 10.0, 20.0, 30.0],
                                            [0.0, 0.0, 0.0, 0.0],
                                            150000.0.*ones(4),
                                            zeros(4),
                                            [50.0.*ones(4) [-2.0, -1.0, 0.0, 1.0] sys.L.*ones(4) -2.5.*ones(4)],
                                            "inlet",
                                            sys.L)
par_g = GasChromatographySimulator.Parameters(sys, prog_g, sub, GasChromatographySimulator.Options())

x = 0.0:par_g.sys.L/1000:par_g.sys.L
grad = Array{Array{Float64,1}}(undef, 4)
gradplot = plot(xlabel="position in m", ylabel="ΔT in °C", legend=:bottomleft)
for i=1:4
    grad[i] = par_g.prog.T_itp.(x, cumsum(par_g.prog.time_steps)[i]).-273.15
    plot!(gradplot, x, grad[i], label="x₀=$(par_g.prog.a_gf[i,2])")
end
gradplot

savefig(gradplot, "grad_x0.svg")

# different gradient profiles - L₀
prog_g = GasChromatographySimulator.Program([0.0, 10.0, 20.0, 30.0],
                                            [0.0, 0.0, 0.0, 0.0],
                                            150000.0.*ones(4),
                                            zeros(4),
                                            [50.0.*ones(4) zeros(4) [9.5, 9.7, 9.9, 10.0] -2.5.*ones(4)],
                                            "inlet",
                                            sys.L)
par_g = GasChromatographySimulator.Parameters(sys, prog_g, sub, GasChromatographySimulator.Options())

x = 0.0:par_g.sys.L/1000:par_g.sys.L
grad = Array{Array{Float64,1}}(undef, 4)
gradplot = plot(xlabel="position in m", ylabel="ΔT in °C", legend=:bottomleft)
for i=1:4
    grad[i] = par_g.prog.T_itp.(x, cumsum(par_g.prog.time_steps)[i]).-273.15
    plot!(gradplot, x, grad[i], label="L₀=$(par_g.prog.a_gf[i,3])")
end
gradplot

savefig(gradplot, "grad_L0.svg")

# different gradient profiles - outlet
prog_g_out = GasChromatographySimulator.Program([0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0],
                                            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                            150000.0.*ones(7),
                                            zeros(7),
                                            [50.0.*ones(7) zeros(7) sys.L.*ones(7) [-10.0, -5.0, -2.5, 0.0, 2.5, 5.0, 10.0]],
                                            "outlet",
                                            sys.L)
par_g_out = GasChromatographySimulator.Parameters(sys, prog_g_out, sub, GasChromatographySimulator.Options(Tcontrol="outlet"))

prog_g_in = GasChromatographySimulator.Program([0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0],
                                            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                            150000.0.*ones(7),
                                            zeros(7),
                                            [50.0.*ones(7) zeros(7) sys.L.*ones(7) [-10.0, -5.0, -2.5, 0.0, 2.5, 5.0, 10.0]],
                                            "inlet",
                                            sys.L)
par_g_in = GasChromatographySimulator.Parameters(sys, prog_g_in, sub, GasChromatographySimulator.Options(Tcontrol="inlet"))


x = 0.0:par_g_in.sys.L/1000:par_g_in.sys.L
T_out = par_g_out.prog.T_itp.(x, 0.0).-273.15
T_in = par_g_in.prog.T_itp.(x, 0.0).-273.15
gradplot = plot(xlabel="position in m", ylabel="ΔT in °C", legend=:bottomleft)
plot!(gradplot, x, T_out, label="Tcontrol=outlet")
plot!(gradplot, x, T_in, label="Tcontrol=inlet")

savefig(gradplot, "grad_Tcontrol.svg")


# self-defined gradient function
L = 10.0
time_steps = [0.0, 60.0, 150.0, 150.0, 120.0]
temp_steps = [40.0, 60.0, 170.0, 300.0, 350.0]
pin_steps = 150000.0.*ones(length(time_steps))
pout_steps = zeros(length(time_steps))
a_gf = [[10.0, 10.0, 30.0, 30.0, 10.0] [1.0, 1.0, 1.0, 2.0, 4.0]]
gradient_function(x) = a_gf[:,1].*sin.(a_gf[:,2].*2*π/L*x)
T_itp = GasChromatographySimulator.temperature_interpolation(time_steps, temp_steps, gradient_function, L)
pin_itp = GasChromatographySimulator.pressure_interpolation(time_steps, pin_steps)
pout_itp = GasChromatographySimulator.pressure_interpolation(time_steps, pout_steps)
prog = GasChromatographySimulator.Program(time_steps, temp_steps, pin_steps, pout_steps, gradient_function, a_gf, T_itp, pin_itp, pout_itp)

x = 0.0:L/1000:L
singradplot = plot(xlabel="position in m", ylabel="ΔT in °C", legend=:bottomleft)
for i=1:5
    grad = prog.T_itp.(x, cumsum(prog.time_steps)[i]).-273.15
    plot!(singradplot, x, grad, label="t[$(i)]")
end
singradplot
# add plot at time inbetween t3 and t4
t = 250.0
plot!(singradplot, x, prog.T_itp.(x, t).-273.15, label="t=$(t)s", line=:dash)
t = 300.0
plot!(singradplot, x, prog.T_itp.(x, t).-273.15, label="t=$(t)s", line=:dash)

savefig(singradplot, "grad_sin.svg")

