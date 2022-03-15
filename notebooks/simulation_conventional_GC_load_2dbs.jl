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
		Pkg.PackageSpec(name="UrlDownload", version="1"),
	Pkg.PackageSpec(url="https://github.com/JanLeppert/GasChromatographySimulator.jl", rev="main"),
Pkg.PackageSpec(url="https://github.com/JanLeppert/GasChromatographyTools.jl", rev="main")		
    ])

    using Plots, PlutoUI, UrlDownload, GasChromatographySimulator, GasChromatographyTools
	md"""
	Packages, simulation\_conventional\_GC\_load\_2dbs.jl, v0.1.0
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
$(Resource("https://raw.githubusercontent.com/JanLeppert/GasChromatographySimulator.jl/main/docs/src/assets/logo.svg"))
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
if typeof(db_file_1) == Nothing
	md"""
	#### !!!Attention!!!
	**Select a 1st database.**
	"""
else
	db_1 = DataFrame(CSV.File(db_file_1["data"]))
	md"""
	1st File: $(db_file_1["name"])
	$(embed_display(db_1))
	"""
end

# ╔═╡ 1b827021-05a6-4fc8-bdce-31302f4a5d0f
if typeof(db_file_2) == Nothing
	md"""
	#### !!!Attention!!!
	**Select a 2nd database.**
	"""
else
	db_2 = DataFrame(CSV.File(db_file_2["data"]))
	md"""
	2nd File: $(db_file_2["name"])
	$(embed_display(db_2))
	"""
end

# ╔═╡ 59e7dff6-7db0-4262-b4f5-ffe8f1a28542
if (@isdefined db_1) && (@isdefined db_2)
	sp_1 = unique(db_1.Phase)
	sp_2 = unique(db_2.Phase)
	sp = GasChromatographyTools.common(sp_1, sp_2)
	common_db = innerjoin(db_1[!, [:Name, :Phase, :Tchar, :thetachar, :DeltaCp, :phi0]], db_2[!, [:Name, :Phase, :Tchar, :thetachar, :DeltaCp, :phi0]], on=[:Name, :Phase], makeunique=true)
	rename!(common_db, [:Name, :Phase, :Tchar_1, :θchar_1, :ΔCp_1, :φ₀_1, :Tchar_2, :θchar_2, :ΔCp_2, :φ₀_2])
	md"""
	Show the common (name/stationary phases) database:
	$(embed_display(common_db))
	"""
elseif @isdefined db_1
	sp = unique(db_1.Phase)
	md"""
	Only data from the 1st database is used.
	"""
elseif @isdefined db_2
	online_db = DataFrame(urldownload("https://github.com/JanLeppert/GasChromatographySimulator.jl/raw/main/data/Database_test.csv"))
	sp_online = unique(online_db.Phase)
	sp_2 = unique(db_2.Phase)
	sp = GasChromatographyTools.common(sp_online, sp_2)
	common_db = innerjoin(online_db[!, [:Name, :Phase, :Tchar, :thetachar, :DeltaCp, :phi0]], db_2[!, [:Name, :Phase, :Tchar, :thetachar, :DeltaCp, :phi0]], on=[:Name, :Phase], makeunique=true)
	rename!(common_db, [:Name, :Phase, :Tchar_1, :θchar_1, :ΔCp_1, :φ₀_1, :Tchar_2, :θchar_2, :ΔCp_2, :φ₀_2])
	md"""
	Common data from the 2nd database and the database from the Github page is used.
	$(embed_display(common_db))
	"""
else
	online_db = DataFrame(urldownload("https://github.com/JanLeppert/GasChromatographySimulator.jl/raw/main/data/Database_test.csv"))
	sp = unique(online_db.Phase)
	md"""
	If no database is selected as 1st database, the database from the Github page is used.
	$(embed_display(online_db))
	"""
end

# ╔═╡ e0669a58-d5ac-4d01-b079-05412b413dda
@bind col_values confirm(GasChromatographyTools.UI_Column(sp))

# ╔═╡ a7e1f0ee-714e-4b97-8741-d4ab5321d5e0
@bind prog_values confirm(GasChromatographyTools.UI_Program(default=("0 60 600 120", "40 40 300 300", "18 18 98 98", "vacuum")))

# ╔═╡ 603e0868-dda5-4063-872f-ce966b6fa8ad
md"""
### Add measured retention times

$(@bind meas_rt_file FilePicker())
"""

# ╔═╡ 5e17ec20-13a1-4f87-a99c-43029907ab9c
if typeof(meas_rt_file) != Nothing
	meas_rt = DataFrame(CSV.File(meas_rt_file["data"]))
	md"""
	Measured retention data:
	
	$(embed_display(meas_rt))
	"""
end

# ╔═╡ 3c856d47-c6c2-40d3-b547-843f9654f48d
md"""
### Plot of local values

Plot $(@bind yy Select(["z", "t", "T", "τ", "σ", "u"]; default="t")) over $(@bind xx Select(["z", "t"]; default="z"))
"""

# ╔═╡ f7f06be1-c8fa-4eee-953f-0d5ea26fafbf
col = GasChromatographySimulator.Column(col_values[1], col_values[2]*1e-3, col_values[3]*1e-6, col_values[4], col_values[5]);

# ╔═╡ 7a00bb54-553f-47f5-b5db-b40d226f4183
begin
	if (@isdefined db_1) && (@isdefined db_2)
		sub_names = GasChromatographyTools.common(
							GasChromatographySimulator.all_solutes(col.sp, db_1), 								GasChromatographySimulator.all_solutes(col.sp, db_2)
							)
	elseif @isdefined db_1
		sub_names = GasChromatographySimulator.all_solutes(col.sp, db_1)
	elseif @isdefined db_2
		sub_names = GasChromatographySimulator.all_solutes(col.sp, db_2)
	else
		sub_names = GasChromatographySimulator.all_solutes(col.sp, online_db)
	end
	@bind sub_values confirm(GasChromatographyTools.UI_Substance(sub_names))
end

# ╔═╡ 0bb1bc3e-9c23-4fbd-9872-fe2e4a2dbdea
prog = GasChromatographyTools.setting_prog(prog_values, col.L);

# ╔═╡ e3277bb4-301a-4a1e-a838-311832b6d6aa
begin
	if @isdefined db_1
		sub_1 = GasChromatographySimulator.load_solute_database(db_1, col.sp, col.gas, sub_values[1], zeros(length(sub_values[1])), zeros(length(sub_values[1])))
	elseif @isdefined online_db
		sub_1 = GasChromatographySimulator.load_solute_database(online_db, col.sp, col.gas, sub_values[1], zeros(length(sub_values[1])), zeros(length(sub_values[1])))
	end
	if @isdefined db_2
		sub_2 = GasChromatographySimulator.load_solute_database(db_2, col.sp, col.gas, sub_values[1], zeros(length(sub_values[1])), zeros(length(sub_values[1])))
	end
end;

# ╔═╡ 115fa61e-8e82-42b2-8eea-9c7e21d97ea8
opt = GasChromatographySimulator.Options(;abstol=1e-8, reltol=1e-5);

# ╔═╡ 85954bdb-d649-4772-a1cd-0bda5d9917e9
begin
	if @isdefined db_1
		par_1 = GasChromatographySimulator.Parameters(col, prog, sub_1, opt)
	elseif @isdefined online_db
		par_1 = GasChromatographySimulator.Parameters(col, prog, sub_1, opt)	
	end
	if @isdefined db_2
		par_2 = GasChromatographySimulator.Parameters(col, prog, sub_2, opt)
	end
end;

# ╔═╡ fdb39284-201b-432f-bff6-986ddbc49a7d
begin
	gr()
	if (@isdefined db_1) && (@isdefined db_2)
		par = par_1
	elseif (@isdefined db_1) || (@isdefined online_db)
		par = par_1
	elseif @isdefined db_2
		par = par_2
	end
	plot_T = GasChromatographySimulator.plot_temperature(par; selector="T(t)")
	plot_p = GasChromatographySimulator.plot_pressure(par.prog)
	xlabel!(plot_p, "")
	plot_F = GasChromatographySimulator.plot_flow(par)
	l = @layout([a{0.65w} [b; c]])
	p_TpF = plot(plot_T, plot_p, plot_F, layout=l, size=(620,300))
	md"""
	## Plot of the program
	
	$(embed_display(p_TpF))
	"""
end

# ╔═╡ 49faa7ea-0f22-45ca-9ab5-338d0db25564
begin
	if @isdefined db_1
		peaklist_1, solution_1 = GasChromatographySimulator.simulate(par_1)
	elseif @isdefined online_db
		peaklist_1, solution_1 = GasChromatographySimulator.simulate(par_1)
	end
	if @isdefined db_2
		peaklist_2, solution_2 = GasChromatographySimulator.simulate(par_2)
	end
	md"""
	## Simulation
	"""
end

# ╔═╡ 14db2d66-eea6-43b1-9caf-2039709d1ddb
if (@isdefined db_1) && (@isdefined db_2)
	md"""
	### Peaklist
	Comparison of simulation with both databases:
	$(embed_display(GasChromatographyTools.compare_peaklist(peaklist_1, peaklist_2)))
	"""
elseif @isdefined db_1
	md"""
	### Peaklist
	$(embed_display(peaklist_1))
	"""
elseif @isdefined db_2
	md"""
	### Peaklist
	$(embed_display(peaklist_2))
	"""
end

# ╔═╡ 9c1e5bc5-8396-4caf-bb80-fe89fea0a5c5
begin
	if @isdefined meas_rt
		if (@isdefined db_1) && (@isdefined db_2)
			compare_1 = GasChromatographyTools.compare_measurement_simulation(meas_rt, peaklist_1)
			compare_2 = GasChromatographyTools.compare_measurement_simulation(meas_rt, peaklist_2)
			RSS_1 = sum(compare_1.ΔtR[isnan.(compare_1.ΔtR).==false].^2)
			RSS_2 = sum(compare_2.ΔtR[isnan.(compare_2.ΔtR).==false].^2)
			md"""
			### Measurement vs. Simulation
			Compare measured retention times with simulation of 1st database:
			$(embed_display(compare_1))
			
			Compare measured retention times with simulation of 2nd database:
			$(embed_display(compare_2))
		
			The residual sum of squares are:
			- 1st database: $(round(RSS_1; digits=2)) s²
			- 2nd database: $(round(RSS_2; digits=2)) s²
			"""
		elseif @isdefined db_1
			compare_1 = GasChromatographyTools.compare_measurement_simulation(meas_rt, peaklist_1)
			RSS_1 = sum(compare_1.ΔtR[isnan.(compare_1.ΔtR).==false].^2)
			md"""
			### Measurement vs. Simulation
			Compare measured retention times with simulation of 1st database:
			$(embed_display(compare_1))
			
			The residual sum of squares is:
			$(round(RSS_1; digits=2)) s²
			"""
		elseif @isdefined db_2
			compare_2 = GasChromatographyTools.compare_measurement_simulation(meas_rt, peaklist_2)
			RSS_2 = sum(compare_2.ΔtR[isnan.(compare_2.ΔtR).==false].^2)
			md"""
			### Measurement vs. Simulation
			Compare measured retention times with simulation of 2nd database:
			$(embed_display(compare_2))
			
			The residual sum of squares is:
			$(round(RSS_2; digits=2)) s²
			"""
		end
	end
end

# ╔═╡ a2287fe8-5aa2-4259-bf7c-f715cc866243
begin
	plotly()
	if (@isdefined db_1) && (@isdefined db_2)
		pchrom = GasChromatographySimulator.plot_chromatogram(peaklist_1, (0,sum(par_1.prog.time_steps)))[1]
		GasChromatographySimulator.plot_chromatogram!(pchrom, peaklist_2, (0,sum(par_2.prog.time_steps)); mirror=true)
		pchrom[1][1][:label] = "1st database"
		pchrom[1][2][:label] = "2nd database"
		ymax = maximum(pchrom[1][1][:y])*1.1
		ymin = minimum(pchrom[1][2][:y])*1.1
	elseif @isdefined db_1
		pchrom = GasChromatographySimulator.plot_chromatogram(peaklist_1, (0,sum(par_1.prog.time_steps)))[1]
		pchrom[1][1][:label] = "1st database"
		ymax = maximum(pchrom[1][1][:y])*1.1
		ymin = 0.0
	elseif @isdefined db_2
		pchrom = GasChromatographySimulator.plot_chromatogram(peaklist_1, (0,sum(par_1.prog.time_steps)))[1]
		GasChromatographySimulator.plot_chromatogram!(pchrom, peaklist_2, (0,sum(par_2.prog.time_steps)); mirror=true)
		pchrom[1][1][:label] = "online database"
		pchrom[1][2][:label] = "2nd database"
		ymax = maximum(pchrom[1][1][:y])*1.1
		ymin = minimum(pchrom[1][2][:y])*1.1
	else
		pchrom = GasChromatographySimulator.plot_chromatogram(peaklist_1, (0,sum(par_1.prog.time_steps)))[1]
		pchrom[1][1][:label] = "online database"
		ymax = maximum(pchrom[1][1][:y])*1.1
		ymin = 0.0
	end
	# add measured data
	if @isdefined meas_rt
		for i=1:length(meas_rt.RT)
			plot!(pchrom, meas_rt.RT[i].*ones(2), [ymin, ymax], c=:orange, line=:solid, label=meas_rt.Name[i])
			annotate!(meas_rt.RT[i], ymax, Plots.text.(meas_rt.Name[i], 10, :orange))
		end
	end
	# finalize
	plot!(pchrom, legend=true, size=(680,400))
	md"""
	### Chromatogram

	$(embed_display(pchrom))
	"""
end

# ╔═╡ 0740f2e6-bce0-4590-acf1-ad4d7cb7c523
begin
	plotly()
	if (@isdefined db_1) && (@isdefined db_2)
		p_1 = GasChromatographyTools.local_plots(xx, yy, solution_1, par_1)
		plot!(p_1, title="1st database")
		p_2 = GasChromatographyTools.local_plots(xx, yy, solution_2, par_2)
		plot!(p_2, title="2nd database")
		p_local = plot(p_1, p_2)
	elseif @isdefined db_1
		p_local = GasChromatographyTools.local_plots(xx, yy, solution_1, par_1)
		plot!(p_local, title="1st database")
	elseif @isdefined db_2
		p_1 = GasChromatographyTools.local_plots(xx, yy, solution_1, par_1)
		plot!(p_1, title="online database")
		p_2 = GasChromatographyTools.local_plots(xx, yy, solution_2, par_2)
		plot!(p_2, title="2nd database")
		p_local = plot(p_1, p_2)
	else
		p_local = GasChromatographyTools.local_plots(xx, yy, solution_1, par_1)
		plot!(p_local, title="online database")
	end
	md"""
	$(embed_display(p_local))
	"""
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
# ╟─84b0de13-2653-4e68-abf4-d98830809844
# ╟─c7c5da86-57ad-4a35-bbaa-39871ee3205d
# ╟─1b827021-05a6-4fc8-bdce-31302f4a5d0f
# ╟─59e7dff6-7db0-4262-b4f5-ffe8f1a28542
# ╟─e0669a58-d5ac-4d01-b079-05412b413dda
# ╟─a7e1f0ee-714e-4b97-8741-d4ab5321d5e0
# ╟─7a00bb54-553f-47f5-b5db-b40d226f4183
# ╟─603e0868-dda5-4063-872f-ce966b6fa8ad
# ╟─5e17ec20-13a1-4f87-a99c-43029907ab9c
# ╟─fdb39284-201b-432f-bff6-986ddbc49a7d
# ╟─49faa7ea-0f22-45ca-9ab5-338d0db25564
# ╟─14db2d66-eea6-43b1-9caf-2039709d1ddb
# ╟─a2287fe8-5aa2-4259-bf7c-f715cc866243
# ╟─9c1e5bc5-8396-4caf-bb80-fe89fea0a5c5
# ╟─3c856d47-c6c2-40d3-b547-843f9654f48d
# ╟─0740f2e6-bce0-4590-acf1-ad4d7cb7c523
# ╟─f7f06be1-c8fa-4eee-953f-0d5ea26fafbf
# ╟─0bb1bc3e-9c23-4fbd-9872-fe2e4a2dbdea
# ╟─e3277bb4-301a-4a1e-a838-311832b6d6aa
# ╟─115fa61e-8e82-42b2-8eea-9c7e21d97ea8
# ╟─85954bdb-d649-4772-a1cd-0bda5d9917e9
# ╟─95e1ca30-9442-4f39-9af0-34bd202fcc24
