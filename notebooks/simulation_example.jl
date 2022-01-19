### A Pluto.jl notebook ###
# v0.17.5

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 7272450e-73b1-11ec-080d-1d1efd32e836
begin
	import Pkg
    # activate the shared project environment
    Pkg.activate(Base.current_project())
	using Plots
	using GasChromatographySimulator
	using Plots
	using PlutoUI
	TableOfContents()
end

# ╔═╡ 9c54bef9-5b70-4cf7-b110-a2f48f5db066
begin
	plotly()
	html"""
	<style>
	  main {
		max-width: 800px;
	  }
	</style>
	"""
end

# ╔═╡ c9246396-3c01-4a36-bc9c-4ed72fd9e325
md"""
# Gas Chromatography Simulator

An example Simulation of a Gas Chromatography (GC) System with a thermal gradient.
"""

# ╔═╡ 8b3011fd-f3df-4ab0-b611-b943d5f3d470
md"""
## Settings
"""

# ╔═╡ 273dcf96-6de4-4380-a00f-ed119bfa13b7
begin
	solute_db_path = "../data" #"/Users/janleppert/Documents/GitHub/GasChromatographySimulator/data/"
	solute_db = "Database_test.csv"
	db = DataFrame(CSV.File(joinpath(solute_db_path, solute_db), header=1, silencewarnings=true))
	sp = unique(db.Phase)
	md"""
	### Solute Database
	$(embed_display(db))
	"""
end

# ╔═╡ 323a769f-55f9-41dd-b8f1-db7928996a52
md"""
## Plot of temperature program

select temperature plot: $(@bind Tplot Select(["T(x,t)", "T(x)", "T(t)"]; default="T(t)"))
"""

# ╔═╡ 3c856d47-c6c2-40d3-b547-843f9654f48d
md"""
### Plot of local values

Plot $(@bind yy Select(["z", "t", "T", "τ", "σ", "u"]; default="t")) over $(@bind xx Select(["z", "t", "T", "τ", "σ", "u"]; default="z"))
"""

# ╔═╡ 95e1ca30-9442-4f39-9af0-34bd202fcc24
md"""
# End
"""

# ╔═╡ 802e4071-b22b-4411-b589-205292aabc75
# functions
begin

"""
    sys_set_UI(sp)

Construct a combined PlutoUI widget for the settings of the GC system with the selectable stationary phases `sp`.    
"""
function sys_set_UI(sp)
	PlutoUI.combine() do Child
		md"""
		### System settings 
		
		``L`` [m]: $(
			Child(NumberField(0.1:0.1:100.0; default=10.0))
		) ``d`` [mm]: $(
			Child(NumberField(0.01:0.01:1.00; default=0.25))
		) ``d_f`` [µm]: $(
			Child(NumberField(0.01:0.01:1.00; default=0.25))
		) stat. phase: $(
			Child(Select(sp; default="SLB5ms"))
		) Gas: $(
			Child(Select(["He", "H2", "N2"]; default="He"))
		) 
			
		"""
	end
end

"""
    prog_set_UI(sp)

Construct a combined PlutoUI widget for the settings of the program of a GC system.    
"""
function prog_set_UI()
	PlutoUI.combine() do Child
		md"""
		### Program settings 
		_Note: Same number of entrys for every text field._
		
		$(
			Child(TextField((50,1); default="0 60 300 300 120"))
		) time steps [s] 
		
		$(
			Child(TextField((50,1); default="40 40 170 300 300"))
		) temperature steps [°C]
		
		$(
			Child(TextField((50,1); default="0 0 40 60 0"))
		) ``ΔT`` steps [°C]
		
		$(
			Child(TextField((50,1); default="-3 -3 -3 -3 -3"))
		) ``α`` steps

		$(
			Child(TextField((50,1); default="18 18 58 98 98"))
		) ``p_{in}`` steps [kPa(g)]

		$(
			Child(TextField((50,1); default="0 0 0 0 0"))
		)``p_{out}`` steps [kPa(a)]
			
		"""
	end
end

"""
    sub_set_UI(sp)

Construct a combined PlutoUI widget for the settings of the substances separated in the simulated GC system with the selectable substances `subs`.    
"""
function sub_set_UI(sol)
	if length(sol)>10
		select_size = 10
	else
		select_size = length(sol)
	end
	PlutoUI.combine() do Child
		md"""
		### Substance settings 
		
		Select Substances: $(
			Child(MultiSelect(sol; default=sol[1:4], size=select_size))
		) 
		
		Injection time [s]: $(
			Child(NumberField(0.0:0.1:100.0; default=0.0))
		) and Injection width [s]: $(
			Child(NumberField(0.00:0.01:10.0; default=0.0))
		) 
		"""
	end
end

"""
    opt_set_UI(sp)

Construct a combined PlutoUI widget for the settings of the options for the simulation.    
"""
function opt_set_UI()
	PlutoUI.combine() do Child
		md"""
		### Option settings 
		
		abstol: 1e $(
			Child(NumberField(-10:1:-3; default=-8))
		) reltol: 1e $(
			Child(NumberField(-8:1:-2; default=-5))
		) Tcontrol: $(
			Child(Select(["inlet", "outlet"]; default="inlet"))
		)
		"""
	end
end

"""
    chromatogram(t::Array{Float64,1}, tR::Array{Float64,1}, τR::Array{Float64,1})

Calculate the chromatogram as a sum of gaussian peaks over the time `t` for peaks centered at retention times `tR` and with peak width `τR`.    
"""
function chromatogram(t::Array{Float64,1}, tR::Array{Float64,1}, τR::Array{Float64,1})
	g(t,tR,τR) = 1/sqrt(2*π*τR^2)*exp(-(t-tR)^2/(2*τR^2))
	chromatograms = Array{Array{Float64,1}}(undef, length(tR))
	for j=1:length(tR)
		chromatograms[j] = g.(t, tR[j], τR[j])
	end
	return sum(chromatograms)
end

"""
    plot_chromatogram(peaklist)

Plot the chromatogram of the peaks listed in `peaklist``.    
"""
function plot_chromatogram(peaklist)
	tMax = maximum(peaklist.tR)*1.05
	t = 0.0:tMax/10000:tMax
	chrom = chromatogram(collect(t), peaklist.tR, peaklist.τR)
	p_chrom = plot(t, chrom, xlabel="time in s", ylabel="abundance", legend=false)
	return p_chrom, t, chrom
end

function plot_flow(par)
	t = 0.0:sum(par.prog.time_steps)/1000.0:sum(par.prog.time_steps)
	F = Array{Float64}(undef, length(t))
	for i=1:length(t)
		F[i] = GasChromatographySimulator.flow(t[i], par.prog.T_itp, par.prog.pin_itp, par.prog.pout_itp, par.sys.L, par.sys.d, par.sys.gas)
	end
	pflow = plot(t, F.*60e6, xlabel="time in s", ylabel="column flow in mL/min", legend=false, size=(800,500))
	return pflow
end

function plot_pressure(par)
	t = 0.0:sum(par.prog.time_steps)/1000.0:sum(par.prog.time_steps)
	pin = Array{Float64}(undef, length(t))
	pout = Array{Float64}(undef, length(t))
	for i=1:length(t)
		pin[i] = par.prog.pin_itp(t[i])
		pout[i] = par.prog.pout_itp(t[i])
	end
	ppress = plot(t, pin, xlabel="time in s", ylabel="pressure in Pa", label="inlet", size=(800,500))
	plot!(ppress, t, pout, label="outlet")
	return ppress
end

function temperature_plot(par, plot_selector)
	if plot_selector=="T(x)"
		nx = 0.0:par.sys.L/1000:par.sys.L
		Tplot = plot(xlabel="x in m", ylabel="T in °C", legend=:top, size=(800,500))
		for i=1:length(par.prog.time_steps)
			T = par.prog.T_itp.(nx, cumsum(par.prog.time_steps)[i]).-273.15
			plot!(Tplot, nx, T, label="t=$(cumsum(par.prog.time_steps)[i])s")
		end
	elseif plot_selector=="T(t)"
		nt = 0.0:sum(par.prog.time_steps)/1000:sum(par.prog.time_steps)
		T0 = par.prog.T_itp.(0.0, nt).-273.15
		Tplot = plot(nt, T0, xlabel="t in s", ylabel="T in °C", legend=:top, label="inlet", c=:red, size=(800,500))
		TL = par.prog.T_itp.(par.sys.L, nt).-273.15
		plot!(Tplot, nt, TL, label="outlet", c=:blue)
	elseif plot_selector=="T(x,t)"
		nx = 0.0:par.sys.L/1000:par.sys.L
		nt = 0.0:sum(par.prog.time_steps)/1000:sum(par.prog.time_steps)
		T = Array{Float64}(undef, length(nx), length(nt))
		for j=1:length(nt)
			for i=1:length(nx)
				T[i,j] = par.prog.T_itp(nx[i], nt[j])-273.15
			end
		end
		Tplot = plot(nx, nt, T', st=:surface, xlabel="x in m", ylabel="t in s", zlabel="T in °C", size=(800,500))
	end
	return Tplot
end

	md"""
	Definition of functions.
	"""
end

# ╔═╡ e0669a58-d5ac-4d01-b079-05412b413dda
@bind sys_values confirm(sys_set_UI(sp))

# ╔═╡ a7e1f0ee-714e-4b97-8741-d4ab5321d5e0
@bind prog_values confirm(prog_set_UI())

# ╔═╡ 3e053ac1-db7b-47c1-b52c-00e26b59912f
@bind opt_values confirm(opt_set_UI())

# ╔═╡ f7f06be1-c8fa-4eee-953f-0d5ea26fafbf
sys = GasChromatographySimulator.System(sys_values[1], sys_values[2]*1e-3, sys_values[3]*1e-6, sys_values[4], sys_values[5]);

# ╔═╡ 7a00bb54-553f-47f5-b5db-b40d226f4183
@bind sub_values confirm(sub_set_UI(GasChromatographySimulator.all_solutes(sys.sp, db)))

# ╔═╡ e3277bb4-301a-4a1e-a838-311832b6d6aa
sub = GasChromatographySimulator.load_solute_database(db, sys.sp, sys.gas, sub_values[1], sub_values[2].*ones(length(sub_values[1])), sub_values[3].*ones(length(sub_values[1])));

# ╔═╡ 115fa61e-8e82-42b2-8eea-9c7e21d97ea8
opt = GasChromatographySimulator.Options(;abstol=10.0^opt_values[1], reltol=10.0^opt_values[2], Tcontrol=opt_values[3]);

# ╔═╡ ee267b33-4086-4e04-9f39-b7f53f2ec920
prog = GasChromatographySimulator.Program(parse.(Float64, split(prog_values[1])),
										parse.(Float64, split(prog_values[2])),
										parse.(Float64, split(prog_values[5])).*1000.0.+101300.0,
										parse.(Float64, split(prog_values[6])).*1000.0,
										parse.(Float64, split(prog_values[3])),
										zeros(length(split(prog_values[1]))),
										sys.L.*ones(length(split(prog_values[1]))),
										parse.(Float64, split(prog_values[4])),
										opt.Tcontrol,
										sys.L
);

# ╔═╡ 85954bdb-d649-4772-a1cd-0bda5d9917e9
par = GasChromatographySimulator.Parameters(sys, prog, sub, opt);

# ╔═╡ 96580972-6ef1-4a88-b872-2f74cba4dbf4
begin
	if Tplot=="T(x,t)"
		md"""
		**_Temperature T(x,t)_**
		
		$(embed_display(temperature_plot(par, Tplot)))
		"""
	elseif Tplot=="T(x)"
		md"""
		**_Temperature T(x)_**
		
		$(embed_display(temperature_plot(par, Tplot)))
		"""
	elseif Tplot=="T(t)"
		md"""
		**_Temperature T(t)_**
		$(embed_display(temperature_plot(par, Tplot)))
		"""
	end
end

# ╔═╡ 1554aee2-41dc-4a39-b6fc-5e02ad4c24e2
md"""
## Plot of pressure program

$(embed_display(plot_pressure(par)))
"""

# ╔═╡ 834a26d2-8f7b-4a00-843f-19e13dc686f2
md"""
## Plot of column flow

$(embed_display(plot_flow(par)))
"""

# ╔═╡ 49faa7ea-0f22-45ca-9ab5-338d0db25564
begin	
	peaklist, solution = GasChromatographySimulator.simulate(par)
	md"""
	## Simulation
	"""
end

# ╔═╡ 14db2d66-eea6-43b1-9caf-2039709d1ddb
md"""
### Peaklist
$(embed_display(peaklist))
"""

# ╔═╡ a2287fe8-5aa2-4259-bf7c-f715cc866243
begin
	pchrom = plot_chromatogram(peaklist)[1]
	md"""
	### Chromatogram

	$(embed_display(pchrom))
	"""
end

# ╔═╡ 48f91bc4-35ce-470a-9b6d-eb3c08da27dc
begin
	
function add_plots(xx, yy, sol, par)
	n = size(sol)[1]

	df_sol = GasChromatographySimulator.sol_extraction(solution, par)
	xvalues = Array{Array{Float64,1}}(undef, n)
	yvalues = Array{Array{Float64,1}}(undef, n)
	
	p_add = plot(legend=false)
	for i=1:n
		if xx=="z"
			xvalues[i] = df_sol.z[i]
			xlabel = "position z in m"
		elseif xx=="t"
			xvalues[i] = df_sol.t[i]
			xlabel = "time t in s"
		elseif xx=="T"
			xvalues[i] = par.prog.T_itp.(df_sol.z[i], df_sol.t[i]).-273.15
			xlabel = "temperature T in °C"
		elseif xx=="τ"
			xvalues[i] = sqrt.(df_sol.τ²[i])
			xlabel = "peak width τ in s"
		elseif xx=="σ"
			xvalues[i] = velocity(df_sol, i, par).*sqrt.(df_sol.τ²[i])
			xlabel = "band width in m"
		elseif xx=="u"
			xvalues[i] = velocity(df_sol, i, par)
			xlabel = "solute velocity in m/s"
		end
		if yy=="z"
			yvalues[i] = df_sol.z[i]
			ylabel = "position z in m"
		elseif yy=="t"
			yvalues[i] = df_sol.t[i]
			ylabel = "time t in s"
		elseif yy=="T"
			yvalues[i] = par.prog.T_itp.(df_sol.z[i], df_sol.t[i]).-273.15
			ylabel = "temperature T in °C"
		elseif yy=="τ"
			yvalues[i] = sqrt.(df_sol.τ²[i])
			ylabel = "peak width in s"
		elseif yy=="σ"
			yvalues[i] = velocity(df_sol, i, par).*sqrt.(df_sol.τ²[i])
			ylabel = "band width in m"
		elseif yy=="u"
			yvalues[i] = velocity(df_sol, i, par)
			ylabel = "solute velocity in m/s"
		end
		plot!(p_add, xvalues[i], yvalues[i], xlabel=xlabel, ylabel=ylabel, label=par.sub[i].name, m=:o)
	end
	return p_add
end

function velocity(df_sol, i, par)
	x = df_sol.z[i]
	t = df_sol.t[i]
	T_itp = par.prog.T_itp
	pin_itp = par.prog.pin_itp
	pout_itp = par.prog.pout_itp
	L = par.sys.L
	d = par.sys.d
	df = par.sys.df
	gas = par.sys.gas
	ΔCp = par.sub[i].ΔCp
	Tchar = par.sub[i].Tchar
	θchar = par.sub[i].θchar
	φ₀ = par.sub[i].φ₀
	u = Array{Float64}(undef, length(x))
	for j=1:length(x)
		u[j] = 1/GasChromatographySimulator.residency(x[j], t[j], T_itp, pin_itp, pout_itp, L, d, df, gas, ΔCp, Tchar, θchar, φ₀)
	end
	return u
end
	md"""
	Definition of more functions.
	"""
end

# ╔═╡ 0740f2e6-bce0-4590-acf1-ad4d7cb7c523
md"""
$(embed_display(add_plots(xx, yy, solution, par)))
"""

# ╔═╡ Cell order:
# ╟─7272450e-73b1-11ec-080d-1d1efd32e836
# ╟─9c54bef9-5b70-4cf7-b110-a2f48f5db066
# ╟─c9246396-3c01-4a36-bc9c-4ed72fd9e325
# ╟─8b3011fd-f3df-4ab0-b611-b943d5f3d470
# ╟─273dcf96-6de4-4380-a00f-ed119bfa13b7
# ╟─e0669a58-d5ac-4d01-b079-05412b413dda
# ╟─a7e1f0ee-714e-4b97-8741-d4ab5321d5e0
# ╟─7a00bb54-553f-47f5-b5db-b40d226f4183
# ╟─3e053ac1-db7b-47c1-b52c-00e26b59912f
# ╠═323a769f-55f9-41dd-b8f1-db7928996a52
# ╟─96580972-6ef1-4a88-b872-2f74cba4dbf4
# ╠═1554aee2-41dc-4a39-b6fc-5e02ad4c24e2
# ╠═834a26d2-8f7b-4a00-843f-19e13dc686f2
# ╠═49faa7ea-0f22-45ca-9ab5-338d0db25564
# ╟─14db2d66-eea6-43b1-9caf-2039709d1ddb
# ╟─a2287fe8-5aa2-4259-bf7c-f715cc866243
# ╟─3c856d47-c6c2-40d3-b547-843f9654f48d
# ╟─0740f2e6-bce0-4590-acf1-ad4d7cb7c523
# ╟─95e1ca30-9442-4f39-9af0-34bd202fcc24
# ╠═802e4071-b22b-4411-b589-205292aabc75
# ╟─48f91bc4-35ce-470a-9b6d-eb3c08da27dc
# ╟─f7f06be1-c8fa-4eee-953f-0d5ea26fafbf
# ╟─ee267b33-4086-4e04-9f39-b7f53f2ec920
# ╟─e3277bb4-301a-4a1e-a838-311832b6d6aa
# ╟─115fa61e-8e82-42b2-8eea-9c7e21d97ea8
# ╟─85954bdb-d649-4772-a1cd-0bda5d9917e9
