### A Pluto.jl notebook ###
# v0.19.26

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
	# online version
	import Pkg
	version = "0.4.3"
	Pkg.activate(mktempdir())
	Pkg.add([
		Pkg.PackageSpec(name="CSV"),
		Pkg.PackageSpec(name="DataFrames"),
		Pkg.PackageSpec(name="GasChromatographySimulator", version=version),
		#Pkg.PackageSpec(name="GasChromatographySimulator", rev="fix_load_solute_database"),
        Pkg.PackageSpec(name="HypertextLiteral"),
		Pkg.PackageSpec(name="Plots"),
		Pkg.PackageSpec(name="PlutoUI"),
		Pkg.PackageSpec(name="UrlDownload"),
    ])
    using CSV, DataFrames,  GasChromatographySimulator, HypertextLiteral, Plots, PlutoUI, UrlDownload
	md"""
	online, Packages, simulation\_example.jl, for GasChromatographySimulator v$(version)
	"""
	
	# local version (database is still downloaded from github)
#=
	import Pkg
	# activate the shared project environment
	Pkg.activate(Base.current_project())
	using CSV, DataFrames,  GasChromatographySimulator, HypertextLiteral, Plots, PlutoUI, UrlDownload
	md"""
	local, Packages, simulation\_example.jl, for GasChromatographySimulator v0.4.1
	"""
=#
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
# $(Resource("https://raw.githubusercontent.com/JanLeppert/GasChromatographySimulator.jl/main/docs/src/assets/logo_b.svg"))
An example Simulation of a Gas Chromatography (GC) System with a thermal gradient.
"""

# ╔═╡ 8b3011fd-f3df-4ab0-b611-b943d5f3d470
md"""
# Settings
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
		db = DataFrame(urldownload("https://raw.githubusercontent.com/JanLeppert/RetentionData/main/Databases/GCSim_database_nonflag.csv"))
	else
		db = DataFrame(CSV.File(db_file["data"], silencewarnings=true, stringtype=String))
	end
	insertcols!(db, 1, :No => collect(1:length(db.Name)))
	sp = unique(db.Phase)
	md"""
	$(embed_display(db))
	"""
end

# ╔═╡ e0669a58-d5ac-4d01-b079-05412b413dda
@bind col_values confirm(GasChromatographySimulator.UI_Column(sp))

# ╔═╡ 560a8d83-0337-45ed-9611-11891af77a82
md"""### Substance category"""

# ╔═╡ 3b453f46-1003-4a58-9ea7-192b00695959
begin
	cat_filter=filter([:Phase]=>(x)-> x.== col_values[4], db)
	
	@bind cat_values confirm(MultiSelect(["all categories"; unique(skipmissing([cat_filter.Cat_1 cat_filter.Cat_2 cat_filter.Cat_3]))]; default=["all categories"]))
end	

# ╔═╡ 323a769f-55f9-41dd-b8f1-db7928996a52
md"""
# Plot of the program

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

# ╔═╡ f7f06be1-c8fa-4eee-953f-0d5ea26fafbf
col = GasChromatographySimulator.Column(col_values[1], col_values[2]*1e-3, col_values[3]*1e-6, col_values[4], col_values[5]);

# ╔═╡ ee5111de-f1c9-452d-a52c-16d7ccca7ba8
begin 	
	if cat_values == ["all categories"]
		@bind sub_values confirm(GasChromatographySimulator.UI_Substance(GasChromatographySimulator.all_solutes(col.sp, db; id=true); default=(1:5, 0.0, 0.0)))
	else	

		dbfilter=
			try 
				filter([:Cat_1]=>(x)-> occursin(string(x), string(cat_values)), db) 
			catch
				try 
					filter([:Cat_2]=>(x)-> occursin(string(x), string(cat_values)), db)
				catch
					try 
						filter([:Cat_3]=>(x)-> occursin(string(x), string(cat_values)), db)
					catch
					end
				end
			end
		@bind sub_values confirm(GasChromatographySimulator.UI_Substance(GasChromatographySimulator.all_solutes(col.sp, dbfilter; id=true); default=(1:1, 0.0, 0.0)))
	end 
end

# ╔═╡ e3277bb4-301a-4a1e-a838-311832b6d6aa
sub = GasChromatographySimulator.load_solute_database(db, col.sp, col.gas, GasChromatographySimulator.pares_No_from_sub_values(sub_values[1]), sub_values[2].*ones(length(sub_values[1])), sub_values[3].*ones(length(sub_values[1])));

# ╔═╡ 3eccee4a-ee64-4910-bddf-faa124d51e1f
function UI_Program(opt; default=("0 60 300 300 120", "40 40 170 300 300", "0 0 40 60 0", "-3 -3 -3 -3 -3", "18 18 58 98 98", "0 0 0 0 0"))
	if opt.control == "Pressure"
		Fpin_text = "p_in steps [kPa(g)]"
	elseif opt.control == "Flow"
		Fpin_text = "Flow steps [mL/min]"
	end
	PlutoUI.combine() do Child
		@htl("""
		<h3>Program settings</h3> 
		<em>Note: Same number of entrys for every text field.</em>
		<ul>
		$(
			Child(TextField((50,1); default=default[1]))
		) time steps [s] 
		
		$(
			Child(TextField((50,1); default=default[2]))
		) temperature steps [°C]
		
		$(
			Child(TextField((50,1); default=default[3]))
		) ΔT steps [°C]
		
		$(
			Child(TextField((50,1); default=default[4]))
		) α steps

		$(
			Child(TextField((50,1); default=default[5]))
		) $(Fpin_text)

		$(
			Child(TextField((50,1); default=default[6]))
		) p_out steps [kPa(a)]
		</ul>	
		""")
	end
end

# ╔═╡ 85e52a8b-2709-4bff-b46b-09a2ea85c2cf
function UI_Options()
	PlutoUI.combine() do Child
		@htl("""
		<h3>Option settings</h3>
		
		abstol: 1e $(
			Child(NumberField(-10:1:-3; default=-8))
		) reltol: 1e $(
			Child(NumberField(-8:1:-2; default=-5))
		) Tcontrol: $(
			Child(Select(["inlet", "outlet"]; default="inlet"))
		)
		<br>
		viscosity model: $(
			Child(Select(["HP", "Blumberg"]; default="Blumberg"))
		) control mode: $(
			Child(Select(["Pressure", "Flow"]; default="Pressure"))
		)
		<br>
		""")
	end
end

# ╔═╡ 3e053ac1-db7b-47c1-b52c-00e26b59912f
@bind opt_values confirm(UI_Options())

# ╔═╡ 115fa61e-8e82-42b2-8eea-9c7e21d97ea8
opt = GasChromatographySimulator.Options(;abstol=10.0^opt_values[1], reltol=10.0^opt_values[2], Tcontrol=opt_values[3], vis=opt_values[4], control=opt_values[5]);

# ╔═╡ a7e1f0ee-714e-4b97-8741-d4ab5321d5e0
@bind prog_values confirm(UI_Program(opt))

# ╔═╡ ee267b33-4086-4e04-9f39-b7f53f2ec920
begin
	if opt.control=="Pressure"
		a = 101300.0
		b = 1000.0
	elseif opt.control=="Flow"
		a = 0.0
		b = 1/(6e7)
	end
	prog = GasChromatographySimulator.Program(parse.(Float64, split(prog_values[1])),
											parse.(Float64, split(prog_values[2])),
											parse.(Float64, split(prog_values[5])).*b.+a,
											parse.(Float64,split(prog_values[6])).*1000.0,
											parse.(Float64, split(prog_values[3])),
											zeros(length(split(prog_values[1]))),
											col.L.*ones(length(split(prog_values[1]))),
											parse.(Float64, split(prog_values[4])),
											opt.Tcontrol,
											col.L
	)
end;

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

# ╔═╡ 629b46e6-0c46-4372-a22a-f3749f6f2a04
par

# ╔═╡ 49faa7ea-0f22-45ca-9ab5-338d0db25564
begin	
	peaklist, solution = GasChromatographySimulator.simulate(par)
	md"""
	# Simulation
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

# ╔═╡ 03ff111e-674a-44c9-89fd-d8f4e255ee53
function export_str(opt_values, col_values, prog_values, pl)
	opt_str_array = ["abstol = 1e$(opt_values[1])", "reltol = 1e$(opt_values[2])", "Tcontrol = $(opt_values[3])", "viscosity = $(opt_values[4])", "control = $(opt_values[5])"]
	opt_str = string(join(opt_str_array, ", "), "\n")
	
	col_str_array = ["L = $(col_values[1]) m", "d = $(col_values[2]) mm", "df = $(col_values[3]) µm", col_values[4], "gas = $(col_values[5])"]
	col_str = string(join(col_str_array, ", "), "\n")

	if opt.control == "Pressure"
		prog_str_array = ["time steps in s: $(prog_values[1])", "temperature steps in °C = $(prog_values[2])", "gradient ΔT in °C = $(prog_values[3])", "gradient α = $(prog_values[4])", "pin in kPa(g) = $(prog_values[5])", "pout in kPa(a) = $(prog_values[6])"]
	elseif opt.control == "Flow"
		prog_str_array = ["Program: $(prog_values[1])", "F = $(prog_values[2]) mL/min", "outlet = $(prog_values[3])"]
	end
	prog_str = string(join(prog_str_array, ", "), "\n")

	header = string(join(names(pl), ", "), "\n")

	pl_array = Array{String}(undef, length(pl.Name))
	for i=1:length(pl.Name)
		pl_array[i] = string(join(Matrix(pl)[i,:], ", "), "\n")
	end
	pl_str = join(pl_array)
	
	export_str = string("Option settings: \n", opt_str, "Column settings: \n", col_str, "Program settings: \n", prog_str, "Peaklist: \n", header, pl_str)
	return export_str
end

# ╔═╡ bef25d30-f654-425b-aad1-926ffff22679
begin
	export_str_ = export_str(opt_values, col_values, prog_values, peaklist)
	md"""
	## Export Results
	Filename: $(@bind result_filename TextField((20,1); default="Result.txt"))
	"""
end

# ╔═╡ 10846818-ea29-4c8a-aec0-20eca782e509
md"""
$(DownloadButton(export_str_, result_filename))
"""

# ╔═╡ Cell order:
# ╠═115b320f-be42-4116-a40a-9cf1b55d39b5
# ╟─9c54bef9-5b70-4cf7-b110-a2f48f5db066
# ╟─c9246396-3c01-4a36-bc9c-4ed72fd9e325
# ╟─8b3011fd-f3df-4ab0-b611-b943d5f3d470
# ╟─7fdd3112-c8df-431c-a514-f386423bca17
# ╟─77a34d59-bc6e-4c5a-ad51-f69903449cb0
# ╠═273dcf96-6de4-4380-a00f-ed119bfa13b7
# ╟─3e053ac1-db7b-47c1-b52c-00e26b59912f
# ╟─e0669a58-d5ac-4d01-b079-05412b413dda
# ╟─a7e1f0ee-714e-4b97-8741-d4ab5321d5e0
# ╟─560a8d83-0337-45ed-9611-11891af77a82
# ╟─3b453f46-1003-4a58-9ea7-192b00695959
# ╟─ee5111de-f1c9-452d-a52c-16d7ccca7ba8
# ╟─323a769f-55f9-41dd-b8f1-db7928996a52
# ╟─fdb39284-201b-432f-bff6-986ddbc49a7d
# ╟─629b46e6-0c46-4372-a22a-f3749f6f2a04
# ╟─49faa7ea-0f22-45ca-9ab5-338d0db25564
# ╟─14db2d66-eea6-43b1-9caf-2039709d1ddb
# ╟─a2287fe8-5aa2-4259-bf7c-f715cc866243
# ╟─3c856d47-c6c2-40d3-b547-843f9654f48d
# ╟─0740f2e6-bce0-4590-acf1-ad4d7cb7c523
# ╟─bef25d30-f654-425b-aad1-926ffff22679
# ╟─10846818-ea29-4c8a-aec0-20eca782e509
# ╟─95e1ca30-9442-4f39-9af0-34bd202fcc24
# ╟─115fa61e-8e82-42b2-8eea-9c7e21d97ea8
# ╟─f7f06be1-c8fa-4eee-953f-0d5ea26fafbf
# ╟─ee267b33-4086-4e04-9f39-b7f53f2ec920
# ╟─e3277bb4-301a-4a1e-a838-311832b6d6aa
# ╟─85954bdb-d649-4772-a1cd-0bda5d9917e9
# ╟─3eccee4a-ee64-4910-bddf-faa124d51e1f
# ╟─85e52a8b-2709-4bff-b46b-09a2ea85c2cf
# ╟─03ff111e-674a-44c9-89fd-d8f4e255ee53
