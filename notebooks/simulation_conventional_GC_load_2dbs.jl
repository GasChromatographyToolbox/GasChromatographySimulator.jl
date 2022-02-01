### A Pluto.jl notebook ###
# v0.17.7

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
	Pkg.PackageSpec(url="https://github.com/JanLeppert/GasChromatographySimulator.jl"),
Pkg.PackageSpec(url="https://github.com/JanLeppert/GasChromatographyTools.jl")		
    ])

    using Plots, PlutoUI, GasChromatographySimulator, GasChromatographyTools
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
# Gas Chromatography Simulator - load two databases

A Simulation of a conventional Gas Chromatography (GC) System (without a thermal gradient).
"""

# ╔═╡ 8b3011fd-f3df-4ab0-b611-b943d5f3d470
md"""
## Settings
"""

# ╔═╡ 84b0de13-2653-4e68-abf4-d98830809844
md"""
### Select Solute Database
first database: $(@bind db_file_1 FilePicker())

second database: $(@bind db_file_2 FilePicker())
"""

# ╔═╡ c7c5da86-57ad-4a35-bbaa-39871ee3205d
begin
	db_1 = DataFrame(CSV.File(db_file_1["data"]))
	sp_1 = unique(db_1.Phase)
	db_2 = DataFrame(CSV.File(db_file_2["data"]))
	sp_2 = unique(db_2.Phase)
	sp = GasChromatographyTools.common(sp_1, sp_2)
	md"""
	1st File: $(db_file_1["name"])
	$(embed_display(db_1))

	2nd File: $(db_file_2["name"])
	$(embed_display(db_2))
	"""
end

# ╔═╡ 59e7dff6-7db0-4262-b4f5-ffe8f1a28542
begin
	common_db = innerjoin(db_1[!, [:Name, :Phase, :Tchar, :thetachar, :DeltaCp, :phi0]], db_2[!, [:Name, :Phase, :Tchar, :thetachar, :DeltaCp, :phi0]], on=[:Name, :Phase], makeunique=true)
	rename!(common_db, [:Name, :Phase, :Tchar_1, :θchar_1, :ΔCp_1, :φ₀_1, :Tchar_2, :θchar_2, :ΔCp_2, :φ₀_2])
	md"""
	Show the common (name/stationary phases) database:
	$(embed_display(common_db))
	"""
end

# ╔═╡ e0669a58-d5ac-4d01-b079-05412b413dda
@bind sys_values confirm(GasChromatographyTools.UI_System(sp))

# ╔═╡ a7e1f0ee-714e-4b97-8741-d4ab5321d5e0
@bind prog_values confirm(GasChromatographyTools.UI_Program_ng())

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

# ╔═╡ f7f06be1-c8fa-4eee-953f-0d5ea26fafbf
sys = GasChromatographySimulator.System(sys_values[1], sys_values[2]*1e-3, sys_values[3]*1e-6, sys_values[4], sys_values[5]);

# ╔═╡ 7a00bb54-553f-47f5-b5db-b40d226f4183
@bind sub_values confirm(GasChromatographyTools.UI_Substance_name(GasChromatographyTools.common(GasChromatographySimulator.all_solutes(sys.sp, db_1), GasChromatographySimulator.all_solutes(sys.sp, db_2))))

# ╔═╡ 0bb1bc3e-9c23-4fbd-9872-fe2e4a2dbdea
prog = GasChromatographyTools.setting_prog(prog_values, sys.L);

# ╔═╡ e3277bb4-301a-4a1e-a838-311832b6d6aa
sub_1 = GasChromatographySimulator.load_solute_database(db_1, sys.sp, sys.gas, sub_values[1], zeros(length(sub_values[1])), zeros(length(sub_values[1])));

# ╔═╡ d201c020-9311-418d-a652-200667d2d64e
sub_2 = GasChromatographySimulator.load_solute_database(db_2, sys.sp, sys.gas, sub_values[1], zeros(length(sub_values[1])), zeros(length(sub_values[1])));

# ╔═╡ 115fa61e-8e82-42b2-8eea-9c7e21d97ea8
opt = GasChromatographySimulator.Options(;abstol=1e-8, reltol=1e-5);

# ╔═╡ 85954bdb-d649-4772-a1cd-0bda5d9917e9
par_1 = GasChromatographySimulator.Parameters(sys, prog, sub_1, opt);

# ╔═╡ fdb39284-201b-432f-bff6-986ddbc49a7d
begin
	gr()
	plot_T = GasChromatographySimulator.plot_temperature(par_1; selector=Tplot)
	if Tplot=="T(x)"
		plot!(plot_T, legend=:bottomleft)
	end
	plot_p = GasChromatographySimulator.plot_pressure(par_1.prog)
	xlabel!(plot_p, "")
	plot_F = GasChromatographySimulator.plot_flow(par_1)
	l = @layout([a{0.65w} [b; c]])
	p_TpF = plot(plot_T, plot_p, plot_F, layout=l, size=(620,300))
	md"""
	$(embed_display(p_TpF))
	"""
end

# ╔═╡ 63d59b94-bea6-4a9a-a0f3-e97602ceeb7e
par_2 = GasChromatographySimulator.Parameters(sys, prog, sub_2, opt);

# ╔═╡ 49faa7ea-0f22-45ca-9ab5-338d0db25564
begin	
	peaklist_1, solution_1 = GasChromatographySimulator.simulate(par_1)
	peaklist_2, solution_2 = GasChromatographySimulator.simulate(par_2)
	md"""
	## Simulation
	"""
end

# ╔═╡ 14db2d66-eea6-43b1-9caf-2039709d1ddb
# new common peaklist name tR1 tR2 DtR tauR1 tauR2 DtauR 
md"""
### Peaklist
$(embed_display(GasChromatographyTools.common_peaklist(peaklist_1, peaklist_2)))
"""

# ╔═╡ a2287fe8-5aa2-4259-bf7c-f715cc866243
begin
	plotly()
	pchrom = GasChromatographySimulator.plot_chromatogram(peaklist_1, (0,sum(par_1.prog.time_steps)))[1]
	GasChromatographySimulator.plot_chromatogram!(pchrom, peaklist_2, (0,sum(par_2.prog.time_steps)); mirror=true)
	md"""
	### Chromatogram

	$(embed_display(pchrom))
	"""
end

# ╔═╡ 0740f2e6-bce0-4590-acf1-ad4d7cb7c523
begin
	plotly()
	p_1 = GasChromatographyTools.local_plots(xx, yy, solution_1, par_1)
	p_2 = GasChromatographyTools.local_plots(xx, yy, solution_2, par_2)
	p_add = plot(p_1, p_2)
end

# ╔═╡ 95e1ca30-9442-4f39-9af0-34bd202fcc24
md"""
# End
"""

# ╔═╡ Cell order:
# ╠═115b320f-be42-4116-a40a-9cf1b55d39b5
# ╟─9c54bef9-5b70-4cf7-b110-a2f48f5db066
# ╟─c9246396-3c01-4a36-bc9c-4ed72fd9e325
# ╟─8b3011fd-f3df-4ab0-b611-b943d5f3d470
# ╟─84b0de13-2653-4e68-abf4-d98830809844
# ╠═c7c5da86-57ad-4a35-bbaa-39871ee3205d
# ╟─59e7dff6-7db0-4262-b4f5-ffe8f1a28542
# ╠═e0669a58-d5ac-4d01-b079-05412b413dda
# ╠═a7e1f0ee-714e-4b97-8741-d4ab5321d5e0
# ╠═7a00bb54-553f-47f5-b5db-b40d226f4183
# ╟─323a769f-55f9-41dd-b8f1-db7928996a52
# ╟─fdb39284-201b-432f-bff6-986ddbc49a7d
# ╟─49faa7ea-0f22-45ca-9ab5-338d0db25564
# ╟─14db2d66-eea6-43b1-9caf-2039709d1ddb
# ╟─a2287fe8-5aa2-4259-bf7c-f715cc866243
# ╟─3c856d47-c6c2-40d3-b547-843f9654f48d
# ╠═0740f2e6-bce0-4590-acf1-ad4d7cb7c523
# ╠═f7f06be1-c8fa-4eee-953f-0d5ea26fafbf
# ╠═0bb1bc3e-9c23-4fbd-9872-fe2e4a2dbdea
# ╠═e3277bb4-301a-4a1e-a838-311832b6d6aa
# ╠═d201c020-9311-418d-a652-200667d2d64e
# ╠═115fa61e-8e82-42b2-8eea-9c7e21d97ea8
# ╠═85954bdb-d649-4772-a1cd-0bda5d9917e9
# ╠═63d59b94-bea6-4a9a-a0f3-e97602ceeb7e
# ╟─95e1ca30-9442-4f39-9af0-34bd202fcc24
