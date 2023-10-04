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

# ╔═╡ 7272450e-73b1-11ec-080d-1d1efd32e836
begin 
	# online version
	import Pkg
	version = "0.4.2"
	Pkg.activate(mktempdir())
	Pkg.add([
		Pkg.PackageSpec(name="CSV"),
		Pkg.PackageSpec(name="DataFrames"),
		Pkg.PackageSpec(name="GasChromatographySimulator", version=version),
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
	local, Packages, simulation\_example\_input\_gradient\_function.jl, for GasChromatographySimulator v0.4.1
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

Test of an arbitrary definitions of gradient functions.
"""

# ╔═╡ 8b3011fd-f3df-4ab0-b611-b943d5f3d470
md"""
# Settings
"""

# ╔═╡ d5332694-7a14-4adc-98e5-1e27deae01ec
md"""
### Solute Database
Load own database: $(@bind own_db CheckBox(default=false))
"""

# ╔═╡ dde6172b-d4e7-45cc-8621-248987c51be6
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

# ╔═╡ 49ac3705-9dbb-48e3-9db3-fc229edb9479
md"""
#### Temperature function $T(x)$
"""

# ╔═╡ 63bd2c94-0626-4eb4-9d2d-368a93a7816d
T(x,a,b) = a.*cos.(b.*π.*x)

# ╔═╡ e5fad38a-1906-4e30-b213-eb68eb3d3265
md"""
!!! note
	Any function of ``x`` with two parameters ``a`` and ``b`` can be defined here.
"""

# ╔═╡ f3d82459-4eca-4612-8fd8-0b1d7b2e400f
md"""### Substance category"""

# ╔═╡ 99ad0726-c4c4-40a7-b729-8ae16a1c13d1
begin
	cat_filter=filter([:Phase]=>(x)-> x.== col_values[4], db)
	
	@bind cat_values confirm(MultiSelect(["all categories"; unique(skipmissing([cat_filter.Cat_1 cat_filter.Cat_2 cat_filter.Cat_3]))]; default=["all categories"]))
end	

# ╔═╡ 3e053ac1-db7b-47c1-b52c-00e26b59912f
@bind opt_values confirm(GasChromatographySimulator.UI_Options())

# ╔═╡ beec4f0f-acfa-41fc-bbf1-67bd5ea9b915
md"""
!!! note
	Default options used for `control="Pressure"` and `vis="Blumberg"`.
"""

# ╔═╡ 323a769f-55f9-41dd-b8f1-db7928996a52
md"""
# Plot of the program

select temperature plot: $(@bind Tplot Select(["T(x,t)", "T(x)", "T(t)"]; default="T(x)"))
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

# ╔═╡ bbebc7a4-6776-4894-be3e-26d8357f26f9
function prog_set_UI_b()
	# Construct a combined PlutoUI widget for the settings of the program of a GC Column. 
	PlutoUI.combine() do Child
		md"""
		### Program settings 
		_Note: Same number of entrys for every text field._
		
		**Gradient function with two parameters ``a`` and ``b``.**
		
		$(
			Child(TextField((50,1); default="0 60 150 150 120"))
		) time steps [s] 
		
		$(
			Child(TextField((50,1); default="40 60 170 300 350"))
		) temperature steps [°C]

		$(
			Child(TextField((50,1); default="50 55 80 110 120"))
		) ``p_{in}`` steps [kPa(g)]

		$(
			Child(TextField((50,1); default="101.3 101.3 101.3 101.3 101.3"))
		)``p_{out}`` steps [kPa(a)]
		
		$(
			Child(TextField((50,1); default="10 10 30 40 20"))
		) ``a`` steps

		$(
			Child(TextField((50,1); default="1 1 2 3 5"))
		) ``b`` steps
			
		"""
	end
end

# ╔═╡ a7e1f0ee-714e-4b97-8741-d4ab5321d5e0
@bind prog_values confirm(prog_set_UI_b())

# ╔═╡ f7f06be1-c8fa-4eee-953f-0d5ea26fafbf
col = GasChromatographySimulator.Column(col_values[1], col_values[2]*1e-3, col_values[3]*1e-6, col_values[4], col_values[5]);

# ╔═╡ 7a00bb54-553f-47f5-b5db-b40d226f4183
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

# ╔═╡ ee267b33-4086-4e04-9f39-b7f53f2ec920
begin
# explicitly call the Program structure to include a self-defined gradient function gf(x)
	time_steps = parse.(Float64, split(prog_values[1]))
	temp_steps = parse.(Float64, split(prog_values[2]))
	pin_steps = parse.(Float64, split(prog_values[3])).*1000.0.+101300.0
	pout_steps = parse.(Float64, split(prog_values[4])).*1000.0
	a_gf = parse.(Float64, split(prog_values[5]))
	b_gf = parse.(Float64, split(prog_values[6]))
	gf(x) = T(x, a_gf, b_gf)
	T_itp = GasChromatographySimulator.temperature_interpolation(time_steps, temp_steps, gf, col.L)
	pin_itp = GasChromatographySimulator.steps_interpolation(time_steps, pin_steps)
    pout_itp = GasChromatographySimulator.steps_interpolation(time_steps, pout_steps)
	
prog = GasChromatographySimulator.Program(time_steps,
										temp_steps,
										pin_steps,
										pout_steps,
										gf,
										[a_gf b_gf],
										T_itp,
										pin_itp,
										pout_itp
)
end;

# ╔═╡ e3277bb4-301a-4a1e-a838-311832b6d6aa
sub = GasChromatographySimulator.load_solute_database(db, col.sp, col.gas, GasChromatographySimulator.pares_No_from_sub_values(sub_values[1]), sub_values[2].*ones(length(sub_values[1])), sub_values[3].*ones(length(sub_values[1])));

# ╔═╡ 115fa61e-8e82-42b2-8eea-9c7e21d97ea8
opt = GasChromatographySimulator.Options(;abstol=10.0^opt_values[1], reltol=10.0^opt_values[2], Tcontrol=opt_values[3]);

# ╔═╡ 85954bdb-d649-4772-a1cd-0bda5d9917e9
par = GasChromatographySimulator.Parameters(col, prog, sub, opt);

# ╔═╡ 96580972-6ef1-4a88-b872-2f74cba4dbf4
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

# ╔═╡ 0c2bb479-31d4-4826-94a8-86e09af9063e
function export_str(opt_values, col_values, prog_values, pl)
	opt_str_array = ["abstol = 1e$(opt_values[1])", "reltol = 1e$(opt_values[2])", "Tcontrol = $(opt_values[3])", "viscosity = Blumberg", "control = Pressure"]
	opt_str = string(join(opt_str_array, ", "), "\n")
	
	col_str_array = ["L = $(col_values[1]) m", "d = $(col_values[2]) mm", "df = $(col_values[3]) µm", col_values[4], "gas = $(col_values[5])"]
	col_str = string(join(col_str_array, ", "), "\n")

	#if opt.control == "Pressure"
		prog_str_array = ["time steps in s: $(prog_values[1])", "temperature steps in °C = $(prog_values[2])", "gradient parameter a = $(prog_values[5])", "gradient parameter b = $(prog_values[6])", "pin in kPa(g) = $(prog_values[3])", "pout in kPa(a) = $(prog_values[4])"]
	#elseif opt.control == "Flow"
	#	prog_str_array = ["Program: $(prog_values[1])", "F = $(prog_values[2]) mL/min", "outlet = $(prog_values[3])"]
	#end
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

# ╔═╡ 3a312086-d89a-484c-a04e-5f572b3acac8
begin
	export_str_ = export_str(opt_values, col_values, prog_values, peaklist)
	md"""
	## Export Results
	Filename: $(@bind result_filename TextField((20,1); default="Result.txt"))
	"""
end

# ╔═╡ 66878953-2c8c-4e95-9e0c-7ac7f82055b3
md"""
$(DownloadButton(export_str_, result_filename))
"""

# ╔═╡ Cell order:
# ╟─7272450e-73b1-11ec-080d-1d1efd32e836
# ╟─9c54bef9-5b70-4cf7-b110-a2f48f5db066
# ╟─c9246396-3c01-4a36-bc9c-4ed72fd9e325
# ╟─8b3011fd-f3df-4ab0-b611-b943d5f3d470
# ╟─d5332694-7a14-4adc-98e5-1e27deae01ec
# ╟─dde6172b-d4e7-45cc-8621-248987c51be6
# ╟─273dcf96-6de4-4380-a00f-ed119bfa13b7
# ╟─e0669a58-d5ac-4d01-b079-05412b413dda
# ╟─a7e1f0ee-714e-4b97-8741-d4ab5321d5e0
# ╟─49ac3705-9dbb-48e3-9db3-fc229edb9479
# ╠═63bd2c94-0626-4eb4-9d2d-368a93a7816d
# ╟─e5fad38a-1906-4e30-b213-eb68eb3d3265
# ╟─f3d82459-4eca-4612-8fd8-0b1d7b2e400f
# ╟─99ad0726-c4c4-40a7-b729-8ae16a1c13d1
# ╟─7a00bb54-553f-47f5-b5db-b40d226f4183
# ╟─3e053ac1-db7b-47c1-b52c-00e26b59912f
# ╟─beec4f0f-acfa-41fc-bbf1-67bd5ea9b915
# ╟─323a769f-55f9-41dd-b8f1-db7928996a52
# ╟─96580972-6ef1-4a88-b872-2f74cba4dbf4
# ╟─49faa7ea-0f22-45ca-9ab5-338d0db25564
# ╟─14db2d66-eea6-43b1-9caf-2039709d1ddb
# ╟─a2287fe8-5aa2-4259-bf7c-f715cc866243
# ╟─3c856d47-c6c2-40d3-b547-843f9654f48d
# ╟─0740f2e6-bce0-4590-acf1-ad4d7cb7c523
# ╟─3a312086-d89a-484c-a04e-5f572b3acac8
# ╟─66878953-2c8c-4e95-9e0c-7ac7f82055b3
# ╟─95e1ca30-9442-4f39-9af0-34bd202fcc24
# ╟─bbebc7a4-6776-4894-be3e-26d8357f26f9
# ╟─f7f06be1-c8fa-4eee-953f-0d5ea26fafbf
# ╟─ee267b33-4086-4e04-9f39-b7f53f2ec920
# ╟─e3277bb4-301a-4a1e-a838-311832b6d6aa
# ╟─115fa61e-8e82-42b2-8eea-9c7e21d97ea8
# ╟─85954bdb-d649-4772-a1cd-0bda5d9917e9
# ╟─0c2bb479-31d4-4826-94a8-86e09af9063e
