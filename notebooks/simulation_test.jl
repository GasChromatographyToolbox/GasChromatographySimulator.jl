### A Pluto.jl notebook ###
# v0.18.1

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
    #import Pkg
    #Pkg.activate(mktempdir())
    #Pkg.add([
    #    Pkg.PackageSpec(name="Plots", version="1"),
    #    Pkg.PackageSpec(name="PlutoUI", version="0.7"),
	#	Pkg.PackageSpec(name="UrlDownload", version="1"),
	#Pkg.PackageSpec(path="/Users/janleppert/Documents/GitHub/GasChromatographySimulator", rev="13-add-const-flow-resp-programmed-flow-mode")
    #])
    using Plots, PlutoUI, UrlDownload, GasChromatographySimulator
	md"""
	Packages, simulation\_test.jl, v0.1.1
	"""
end

# ╔═╡ ed730122-8b5a-40d7-bc0b-f54d88b33fff
pwd()

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
An example Simulation of a Gas Chromatography (GC) System with a thermal gradient.
"""

# ╔═╡ 8b3011fd-f3df-4ab0-b611-b943d5f3d470
md"""
## Settings
"""

# ╔═╡ 7fdd3112-c8df-431c-a514-f386423bca17
md"""
### Solute Database
Load own database: $(@bind own_db CheckBox(default=false))
"""

# ╔═╡ 77a34d59-bc6e-4c5a-ad51-f69903449cb0
if own_db == true
	md"""
	$(@bind db_file FilePicker())
	"""
end

# ╔═╡ 273dcf96-6de4-4380-a00f-ed119bfa13b7
begin
	if own_db == false
		db = DataFrame(urldownload("https://github.com/JanLeppert/GasChromatographySimulator.jl/raw/main/data/Database_test.csv"))
	else
		db = DataFrame(CSV.File(db_file["data"]))
	end
	sp = unique(db.Phase)
	md"""
	$(embed_display(db))
	"""
end

# ╔═╡ e0669a58-d5ac-4d01-b079-05412b413dda
@bind col_values confirm(GasChromatographySimulator.UI_Column(sp))

# ╔═╡ a7e1f0ee-714e-4b97-8741-d4ab5321d5e0
#@bind prog_values confirm(GasChromatographySimulator.UI_Program())

# ╔═╡ 3e053ac1-db7b-47c1-b52c-00e26b59912f
#@bind opt_values confirm(GasChromatographySimulator.UI_Options())

# ╔═╡ 93ff2119-a1c8-4aae-aaa2-84928394dcb0
opt = GasChromatographySimulator.Options(;abstol=1e-8, reltol=1e-5, Tcontrol="inlet", control="Flow", vis="Blumberg")

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
col = GasChromatographySimulator.Column(col_values[1], col_values[2]*1e-3, col_values[3]*1e-6, col_values[4], col_values[5]);

# ╔═╡ 502f7507-7b8a-46dc-8b6b-42dd28c5c4d0
prog = GasChromatographySimulator.Program([0, 60, 300, 300, 120], 
											[40, 40, 170, 300, 300],
											[1, 1, 1, 1, 1]./6e7, 
											zeros(5), 
											[0, 0, 40, 60, 0], 
											zeros(5), 
											col.L.*ones(5), 
											[-3, -3, -3, -3, -3],
											opt.Tcontrol,
											col.L)

# ╔═╡ 7a00bb54-553f-47f5-b5db-b40d226f4183
@bind sub_values confirm(GasChromatographySimulator.UI_Substance(GasChromatographySimulator.all_solutes(col.sp, db); default=(1:5, 0.0, 0.0)))

# ╔═╡ ee267b33-4086-4e04-9f39-b7f53f2ec920
#=prog = GasChromatographySimulator.Program(parse.(Float64, split(prog_values[1])),
										parse.(Float64, split(prog_values[2])),
										parse.(Float64, split(prog_values[5])).*1000.0.+101300.0,
										parse.(Float64, split(prog_values[6])).*1000.0,
										parse.(Float64, split(prog_values[3])),
										zeros(length(split(prog_values[1]))),
										col.L.*ones(length(split(prog_values[1]))),
										parse.(Float64, split(prog_values[4])),
										opt.Tcontrol,
										col.L
);=#

# ╔═╡ e3277bb4-301a-4a1e-a838-311832b6d6aa
sub = GasChromatographySimulator.load_solute_database(db, col.sp, col.gas, sub_values[1], sub_values[2].*ones(length(sub_values[1])), sub_values[3].*ones(length(sub_values[1])));

# ╔═╡ 115fa61e-8e82-42b2-8eea-9c7e21d97ea8
#opt = GasChromatographySimulator.Options(;abstol=10.0^opt_values[1], reltol=10.0^opt_values[2], Tcontrol=opt_values[3]);

# ╔═╡ 85954bdb-d649-4772-a1cd-0bda5d9917e9
par = GasChromatographySimulator.Parameters(col, prog, sub, opt);

# ╔═╡ fdb39284-201b-432f-bff6-986ddbc49a7d
begin
	gr()
	plot_T = GasChromatographySimulator.plot_temperature(par; selector=Tplot)
	if Tplot=="T(x)"
		plot!(plot_T, legend=:bottomleft)
	end
	plot_p = GasChromatographySimulator.plot_pressure(par)
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
	GasChromatographySimulator.local_plots(xx, yy, solution, par)
end

# ╔═╡ 95e1ca30-9442-4f39-9af0-34bd202fcc24
md"""
# End
"""

# ╔═╡ Cell order:
# ╠═115b320f-be42-4116-a40a-9cf1b55d39b5
# ╠═ed730122-8b5a-40d7-bc0b-f54d88b33fff
# ╠═9c54bef9-5b70-4cf7-b110-a2f48f5db066
# ╠═c9246396-3c01-4a36-bc9c-4ed72fd9e325
# ╠═8b3011fd-f3df-4ab0-b611-b943d5f3d470
# ╠═7fdd3112-c8df-431c-a514-f386423bca17
# ╠═77a34d59-bc6e-4c5a-ad51-f69903449cb0
# ╠═273dcf96-6de4-4380-a00f-ed119bfa13b7
# ╠═e0669a58-d5ac-4d01-b079-05412b413dda
# ╠═a7e1f0ee-714e-4b97-8741-d4ab5321d5e0
# ╠═502f7507-7b8a-46dc-8b6b-42dd28c5c4d0
# ╠═7a00bb54-553f-47f5-b5db-b40d226f4183
# ╠═3e053ac1-db7b-47c1-b52c-00e26b59912f
# ╠═93ff2119-a1c8-4aae-aaa2-84928394dcb0
# ╠═323a769f-55f9-41dd-b8f1-db7928996a52
# ╠═fdb39284-201b-432f-bff6-986ddbc49a7d
# ╠═49faa7ea-0f22-45ca-9ab5-338d0db25564
# ╠═14db2d66-eea6-43b1-9caf-2039709d1ddb
# ╠═a2287fe8-5aa2-4259-bf7c-f715cc866243
# ╠═3c856d47-c6c2-40d3-b547-843f9654f48d
# ╠═0740f2e6-bce0-4590-acf1-ad4d7cb7c523
# ╠═f7f06be1-c8fa-4eee-953f-0d5ea26fafbf
# ╠═ee267b33-4086-4e04-9f39-b7f53f2ec920
# ╠═e3277bb4-301a-4a1e-a838-311832b6d6aa
# ╠═115fa61e-8e82-42b2-8eea-9c7e21d97ea8
# ╠═85954bdb-d649-4772-a1cd-0bda5d9917e9
# ╟─95e1ca30-9442-4f39-9af0-34bd202fcc24
