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

# ╔═╡ 115b320f-be42-4116-a40a-9cf1b55d39b5
begin
    import Pkg
    Pkg.activate(mktempdir())
    Pkg.add([
        Pkg.PackageSpec(name="Plots", version="1"),
        Pkg.PackageSpec(name="PlutoUI", version="0.7"),
		Pkg.PackageSpec(url="https://github.com/JanLeppert/GasChromatographySimulator.jl")
    ])
    using Plots, PlutoUI, GasChromatographySimulator
	md"""
	Packages
	"""
end

# ╔═╡ 9c54bef9-5b70-4cf7-b110-a2f48f5db066
begin
	#plotly()
	gr()
	html"""
	<style>
	  main {
		max-width: 800px;
	  }
	</style>
	"""
	TableOfContents()
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
	solute_db_path = "../data"
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
## Plot of the program

select temperature plot: $(@bind Tplot Select(["T(x,t)", "T(x)", "T(t)"]; default="T(t)"))
"""

# ╔═╡ 3c856d47-c6c2-40d3-b547-843f9654f48d
md"""
### Plot of local values

Plot $(@bind yy Select(["z", "t", "T", "τ", "σ", "u"]; default="t")) over $(@bind xx Select(["z", "t", "T", "τ", "σ", "u"]; default="z"))
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



	md"""
	Definition of UI functions.
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

# ╔═╡ fdb39284-201b-432f-bff6-986ddbc49a7d
begin
	gr()
	plot_T = GasChromatographySimulator.plot_temperature(par; selector=Tplot)
	if Tplot=="T(x)"
		plot!(plot_T, legend=:bottomleft)
	end
	plot_p = GasChromatographySimulator.plot_pressure(par.prog)
	xlabel!(plot_p, "")
	plot_F = GasChromatographySimulator.plot_flow(par)
	l = @layout([a{0.65w} [b; c]])
	p_TpF = plot(plot_T, plot_p, plot_F, layout=l)
	md"""
	$(embed_display(p_TpF))
	"""
end

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

# ╔═╡ a2287fe8-5aa2-4259-bf7c-f715cc866243
begin
	plotly()
	pchrom = GasChromatographySimulator.plot_chromatogram(peaklist, (0,sum(par.prog.time_steps)))[1]
	md"""
	### Chromatogram

	$(embed_display(pchrom))
	"""
end

# ╔═╡ 0740f2e6-bce0-4590-acf1-ad4d7cb7c523
begin
	plotly()
	add_plots(xx, yy, solution, par)
end

# ╔═╡ 95e1ca30-9442-4f39-9af0-34bd202fcc24
md"""
# End
"""

# ╔═╡ Cell order:
# ╟─115b320f-be42-4116-a40a-9cf1b55d39b5
# ╟─9c54bef9-5b70-4cf7-b110-a2f48f5db066
# ╟─c9246396-3c01-4a36-bc9c-4ed72fd9e325
# ╟─8b3011fd-f3df-4ab0-b611-b943d5f3d470
# ╟─273dcf96-6de4-4380-a00f-ed119bfa13b7
# ╟─e0669a58-d5ac-4d01-b079-05412b413dda
# ╟─a7e1f0ee-714e-4b97-8741-d4ab5321d5e0
# ╟─7a00bb54-553f-47f5-b5db-b40d226f4183
# ╟─3e053ac1-db7b-47c1-b52c-00e26b59912f
# ╟─323a769f-55f9-41dd-b8f1-db7928996a52
# ╟─fdb39284-201b-432f-bff6-986ddbc49a7d
# ╟─49faa7ea-0f22-45ca-9ab5-338d0db25564
# ╟─14db2d66-eea6-43b1-9caf-2039709d1ddb
# ╟─a2287fe8-5aa2-4259-bf7c-f715cc866243
# ╟─3c856d47-c6c2-40d3-b547-843f9654f48d
# ╟─0740f2e6-bce0-4590-acf1-ad4d7cb7c523
# ╟─802e4071-b22b-4411-b589-205292aabc75
# ╟─48f91bc4-35ce-470a-9b6d-eb3c08da27dc
# ╟─f7f06be1-c8fa-4eee-953f-0d5ea26fafbf
# ╟─ee267b33-4086-4e04-9f39-b7f53f2ec920
# ╟─e3277bb4-301a-4a1e-a838-311832b6d6aa
# ╟─115fa61e-8e82-42b2-8eea-9c7e21d97ea8
# ╟─85954bdb-d649-4772-a1cd-0bda5d9917e9
# ╟─95e1ca30-9442-4f39-9af0-34bd202fcc24
