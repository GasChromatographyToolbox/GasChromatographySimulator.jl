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
	local, Packages, simulation\_conventional\_GC\_TL.jl, for GasChromatographySimulator v0.4.1
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
# 
$(Resource("https://raw.githubusercontent.com/JanLeppert/GasChromatographySimulator.jl/main/docs/src/assets/logo_b.svg"))
A Simulation of a conventional Gas Chromatography (GC) System (without a thermal gradient) consisting of a column with a program and a transfer line with a fixed temperature.
"""

# ╔═╡ 8b3011fd-f3df-4ab0-b611-b943d5f3d470
md"""
# Settings
"""

# ╔═╡ 17966423-96f5-422f-9734-4ab0edab86bd
md"""
### Solute Database
Load own database: $(@bind own_db CheckBox(default=false))
"""

# ╔═╡ 3a076b77-5cd6-4e10-9714-7553d2822806
if own_db == true
	md"""
	$(@bind db_file FilePicker())
	"""
end

# ╔═╡ a0968f4b-b249-4488-a11b-dc109c68150f
begin
	if own_db == false
		db = DataFrame(urldownload("https://raw.githubusercontent.com/JanLeppert/RetentionData/main/Databases/GCSim_database_nonflag.csv"))
	else
		db = DataFrame(CSV.File(db_file["data"], silencewarnings=true))
	end
	insertcols!(db, 1, :No => collect(1:length(db.Name)))
	sp = unique(db.Phase)
	md"""
	$(embed_display(db))
	"""
end

# ╔═╡ e0669a58-d5ac-4d01-b079-05412b413dda
@bind col_values confirm(GasChromatographySimulator.UI_Column(sp))

# ╔═╡ 51eb4859-20b9-4cac-bde4-ef30c6fea59d
md"""
### Program settings

Number of ramps: $(@bind n_ramp confirm(NumberField(0:1:100; default=1)))
"""

# ╔═╡ 6cc59301-fcc6-4bf4-b9ac-5d1b91560e4a
md"""### Substance category"""

# ╔═╡ ae93afbc-051f-47b2-a6cc-4b9f736ef9a8
begin
	cat_filter=filter([:Phase]=>(x)-> x.== col_values[4], db)
	
	@bind cat_values confirm(MultiSelect(["all categories"; unique(skipmissing([cat_filter.Cat_1 cat_filter.Cat_2 cat_filter.Cat_3]))]; default=["all categories"]))
end	

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

# ╔═╡ 7a00bb54-553f-47f5-b5db-b40d226f4183
begin 	
	if cat_values == ["all categories"]
		@bind sub_values confirm(GasChromatographySimulator.UI_Substance(GasChromatographySimulator.all_solutes(col.sp, db; id=true); default=(1:5,)))
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
		@bind sub_values confirm(GasChromatographySimulator.UI_Substance(GasChromatographySimulator.all_solutes(col.sp, dbfilter; id=true); default=(1:1,)))
	end 
end

# ╔═╡ e3277bb4-301a-4a1e-a838-311832b6d6aa
sub = GasChromatographySimulator.load_solute_database(db, col.sp, col.gas, GasChromatographySimulator.pares_No_from_sub_values(sub_values[1]), zeros(length(sub_values[1])), zeros(length(sub_values[1])));

# ╔═╡ 8c831fdb-0bfa-4f36-b720-e82fcf5d2427
function UI_Options()
	PlutoUI.combine() do Child
		@htl("""
		<h3>Option settings</h3>
		
		viscosity model: $(
			Child(Select(["HP", "Blumberg"]; default="Blumberg"))
		) control mode: $(
			Child(Select(["Pressure", "Flow"]; default="Flow"))
		)
		""")
	end
end

# ╔═╡ 0e2eba31-5643-4d10-9ed4-4454ec28df12
@bind opt_values confirm(UI_Options())

# ╔═╡ 115fa61e-8e82-42b2-8eea-9c7e21d97ea8
opt = GasChromatographySimulator.Options(;abstol=1e-8, reltol=1e-5, ng=true, vis=opt_values[1], control=opt_values[2]);

# ╔═╡ c26c1ea7-575f-495d-86bc-987aca991664
function UI_TP(n_ramp)
	PlutoUI.combine() do Child
		@htl("""
		<table>
			<tr>
				<th>ramp [°C/min]</th>
				<th>T [°C]</th>
				<th>hold [min]</th>
			</tr>
			<tr>
				<td></td>
				<td><center>$(Child(NumberField(0.0:0.1:400.0; default=40.0)))</center></td>
				<td><center>$(Child(NumberField(0.0:0.1:100.0; default=1.0)))</center></td>
			</tr>
			$([
				@htl("
					<tr>
						<td><center>$(Child(NumberField(0.0:0.1:100.0; default=5.0)))</center></td>
						<td><center>$(Child(NumberField(0.0:0.1:400.0; default=((i+1)*100.0))))</center></td>
						<td><center>$(Child(NumberField(0.0:0.1:100.0; default=1.0)))</center></td>
					</tr>
				")
				for i=1:n_ramp
			])
		</table>
		""")
	end
end

# ╔═╡ 50f1bd7f-a479-453d-a8ea-57c37a4e330c
function UI_Program(n_ramp, opt)
	if opt.control == "Pressure"
		PlutoUI.combine() do Child
			@htl("""
			<ul>
				temperature program:
				<br>
				$(Child(UI_TP(n_ramp)))
				<br>
				inlet pressure [kPa(g)]:
				$(Child(NumberField(0.0:0.1:1000.0; default=100.0)))
				<br>
				outlet pressure:
				$(Child(Select(["vacuum", "atmosphere"]; default="vacuum")))
			</ul>
			""")
		end
	elseif opt.control == "Flow"
		PlutoUI.combine() do Child
			@htl("""
			<ul>
				temperature program:
				<br>
				$(Child(UI_TP(n_ramp)))
				<br>
				flow [mL/min]:
				$(Child(NumberField(0.0:0.1:100.0; default=1.0)))
				<br>
				outlet pressure:
				$(Child(Select(["vacuum", "atmosphere"]; default="atmosphere")))
			</ul>
			""")
		end
	end
end

# ╔═╡ 83851755-fe6c-4751-aa8e-3226e0fd50da
@bind prog_values confirm(UI_Program(n_ramp, opt))

# ╔═╡ 13af6325-aa2c-44e2-8e2b-d05d8c5ae4a7
"""
    UI_Transferline(sp)

Construct a combined PlutoUI widget for the settings of a transfer line (Column at constant temperature). 
	
**UI fields**
* ``L``: transfer line length in m.
* ``d``: transfer line diameter in mm.
* ``d_f``: film thickness in µm.
* stat. phase: stationary phase of the transfer line
* T: temperature in °C of the transfer line
"""
function UI_Transferline(sp; default=(0.25, 0.25, 0.25, 1, 300.0))
		PlutoUI.combine() do Child
			@htl("""
			<h3>Transfer line settings</h3>
			L [m]: $(
				Child(NumberField(0.1:0.1:100.0; default=default[1]))
			)  d [mm]: $(
				Child(NumberField(0.01:0.01:1.00; default=default[2]))
			)  d_f [µm]: $(
				Child(NumberField(0.01:0.01:1.00; default=default[3]))
			)  stat. phase: $(
				Child(Select([""; sp]; default=sp[default[4]]))
			)  T [°C]: $(
				Child(NumberField(0.0:1.0:500.0; default=default[5]))
			) 
			
			""")
	end
end

# ╔═╡ ced08ac1-8db5-42ae-b8f4-71e154ab2697
#@bind tl_values confirm(GasChromatographySimulator.UI_Transferline(col_values[4]))
@bind tl_values confirm(UI_Transferline(col_values[4]))

# ╔═╡ 331846bd-406f-44e5-94db-78cf0a2594f6
tl = GasChromatographySimulator.Column(tl_values[1], tl_values[2]*1e-3, tl_values[3]*1e-6, tl_values[4], col_values[5]);

# ╔═╡ 5f5e16ec-2730-4a17-bd64-a751426a033f
begin # Program of the column and transferline
	TP_vector = collect(prog_values[1])
	time_steps, temp_steps = GasChromatographySimulator.conventional_program(TP_vector)

	if prog_values[3] == "vacuum"
		pout_steps = zeros(length(time_steps))	
	elseif prog_values[3] == "atmosphere"
		pout_steps = 101300.0.*ones(length(time_steps))
	end

	# for calculation of pressure between column and transferline
	gf(x) = GasChromatographySimulator.gradient(x, [0.0])
	T_itp_col = GasChromatographySimulator.temperature_interpolation(time_steps, temp_steps, gf, col.L)
	κL_col(t) = GasChromatographySimulator.flow_restriction(col.L, t, T_itp_col, col.d, col.gas; ng=opt.ng, vis=opt.vis)
	T_itp_tl = GasChromatographySimulator.temperature_interpolation(time_steps, tl_values[5].*ones(length(time_steps)), gf, tl.L)
	κL_tl(t) = GasChromatographySimulator.flow_restriction(tl.L, t, T_itp_tl, tl.d, tl.gas; ng=opt.ng, vis=opt.vis)
	p_o_col = Array{Float64}(undef, length(time_steps))
	if opt.control=="Pressure"
		a = 101300.0
		b = 1000.0
		for i=1:length(time_steps) # calc pressure between col and tl
			p_o_col[i] = sqrt(((a+b*prog_values[2])^2*κL_tl(time_steps[i]) + pout_steps[1]^2*κL_col(time_steps[i]))/(κL_col(time_steps[i]) + κL_tl(time_steps[i])))
		end
		Fpin_steps = (a + b * prog_values[2]).*ones(length(time_steps))
		Fpin_tl = p_o_col
	elseif opt.control=="Flow"
		a = 0.0
		b = 1/(6e7)
		for i=1:length(time_steps) # calc pressure between col and tl
			p_o_col[i] = sqrt(pout_steps[i]^2 + 256/π * 101300/273.15 * GasChromatographySimulator.viscosity(tl.L, time_steps[i], T_itp_tl, tl.gas; vis=opt.vis) * T_itp_tl(tl.L, time_steps[i]) * tl.L / tl.d^4 * prog_values[2]*b)
		end
		Fpin_steps = (a + b * prog_values[2]).*ones(length(time_steps))
		Fpin_tl = Fpin_steps
	end
	
	# with additional transferline we have to calculate the outlet pressure of the column
	# (p_i,col^2-p_o,col^2)/κ_col - (p_i,tl^2 - p_o,tl^2)/κ_tl = 0
	# 
	
	prog_col = GasChromatographySimulator.Program( 	time_steps,
												temp_steps,
												Fpin_steps,
												p_o_col,
												col.L
												)

	prog_tl = GasChromatographySimulator.Program(time_steps,
												tl_values[5].*ones(length(time_steps)),
												Fpin_tl,
												pout_steps,
												tl.L)
end;

# ╔═╡ 85954bdb-d649-4772-a1cd-0bda5d9917e9
par_col = GasChromatographySimulator.Parameters(col, prog_col, sub, opt);

# ╔═╡ fdb39284-201b-432f-bff6-986ddbc49a7d
begin
	gr()
	plot_T = GasChromatographySimulator.plot_temperature(par_col; selector="T(t)")
	plot_p = GasChromatographySimulator.plot_pressure(par_col)
	xlabel!(plot_p, "")
	plot_F = GasChromatographySimulator.plot_flow(par_col)
	l = @layout([a{0.65w} [b; c]])
	p_TpF = plot(plot_T, plot_p, plot_F, layout=l, size=(620,300))
	md"""
	# Plot of the program

	Pressure/flow for the column is shown, not for the transfer line.
	
	$(embed_display(p_TpF))
	"""
end

# ╔═╡ 49faa7ea-0f22-45ca-9ab5-338d0db25564
begin
	# for definition of par_col and others see cells after section "End"
	peaklist_col, solution_col = GasChromatographySimulator.simulate(par_col)
	sub_tl = GasChromatographySimulator.load_solute_database(db, tl.sp, tl.gas, sub_values[1], peaklist_col.tR, peaklist_col.τR)
	# ATTENTION: handover of the results from the column as initial values to the transferline is not correct -> replace this notebook using GasChromatographySystems.jl
	par_tl = GasChromatographySimulator.Parameters(tl, prog_tl, sub_tl, opt)
	peaklist_tl, solution_tl = GasChromatographySimulator.simulate(par_tl)
	md"""
	# Simulation
	"""
end

# ╔═╡ 14db2d66-eea6-43b1-9caf-2039709d1ddb
md"""
### Peaklist
* After the column:

$(embed_display(peaklist_col))

* After the transfer line:

$(embed_display(peaklist_tl))
"""

# ╔═╡ a2287fe8-5aa2-4259-bf7c-f715cc866243
begin
	plotly()
	pchrom = GasChromatographySimulator.plot_chromatogram(peaklist_tl, (0,sum(par_tl.prog.time_steps)))[1]
	md"""
	### Chromatogram

	$(embed_display(pchrom))
	"""
end

# ╔═╡ ef0f33bd-8e7c-4b88-9566-df8b101828c0
begin
	plotly()
	GasChromatographySimulator.local_plots(xx, yy, solution_tl, par_tl)
	plot!(title="transferline")
end

# ╔═╡ 0740f2e6-bce0-4590-acf1-ad4d7cb7c523
begin
	plotly()
	GasChromatographySimulator.local_plots(xx, yy, solution_col, par_col)
	plot!(title="column")
end

# ╔═╡ 3a54cad5-141f-4258-b767-49b36b6ceccc
function export_str(opt_values, col_values, tl_values, prog_values, pl_col, pl_tl)
	opt_str_array = ["viscosity = $(opt_values[1])", "control = $(opt_values[2])"]
	opt_str = string(join(opt_str_array, ", "), "\n")
	
	col_str_array = ["L = $(col_values[1]) m", "d = $(col_values[2]) mm", "df = $(col_values[3]) µm", col_values[4], "gas = $(col_values[5])"]
	col_str = string(join(col_str_array, ", "), "\n")

	if tl_values[4] == ""
		tl_str_array = ["L = $(tl_values[1]) m", "d = $(tl_values[2]) mm", "df = $(tl_values[3]) µm", "no phase", "T = $(tl_values[5])°C"]
	else
		tl_str_array = ["L = $(tl_values[1]) m", "d = $(tl_values[2]) mm", "df = $(tl_values[3]) µm", tl_values[4], "T = $(tl_values[5])°C"]
	end
	tl_str = string(join(tl_str_array, ", "), "\n")

	if opt.control == "Pressure"
		prog_str_array = ["Program: $(prog_values[1])", "pin = $(prog_values[2]) kPa(g)", "outlet = $(prog_values[3])"]
	elseif opt.control == "Flow"
		prog_str_array = ["Program: $(prog_values[1])", "F = $(prog_values[2]) mL/min", "outlet = $(prog_values[3])"]
	end
	prog_str = string(join(prog_str_array, ", "), "\n")

	header_col = string(join(names(pl_col), ", "), "\n")

	pl_col_array = Array{String}(undef, length(pl_col.Name))
	for i=1:length(pl_col.Name)
		pl_col_array[i] = string(join(Matrix(pl_col)[i,:], ", "), "\n")
	end
	pl_col_str = join(pl_col_array)

	header_tl = string(join(names(pl_tl), ", "), "\n")

	pl_tl_array = Array{String}(undef, length(pl_tl.Name))
	for i=1:length(pl_tl.Name)
		pl_tl_array[i] = string(join(Matrix(pl_tl)[i,:], ", "), "\n")
	end
	pl_tl_str = join(pl_tl_array)
	
	export_str = string("Option settings: \n", opt_str, "Column settings: \n", col_str, "Transfer line settings: \n", tl_str, "Program settings: \n", prog_str, "Peaklist column: \n", header_col, pl_col_str, "Peaklist transfer line: \n", header_tl, pl_tl_str)
	return export_str
end

# ╔═╡ 7a5a4bfa-6640-494e-88e6-1e6d5c75bb77
begin
	export_str_ = export_str(opt_values, col_values, tl_values, prog_values, peaklist_col, peaklist_tl)
	md"""
	## Export Results
	Filename: $(@bind result_filename TextField((20,1); default="Result.txt"))
	"""
end

# ╔═╡ 5211098f-a4b9-4d1e-aaf2-d208b51a2d85
md"""
$(DownloadButton(export_str_, result_filename))
"""

# ╔═╡ Cell order:
# ╠═115b320f-be42-4116-a40a-9cf1b55d39b5
# ╟─9c54bef9-5b70-4cf7-b110-a2f48f5db066
# ╟─c9246396-3c01-4a36-bc9c-4ed72fd9e325
# ╟─8b3011fd-f3df-4ab0-b611-b943d5f3d470
# ╟─17966423-96f5-422f-9734-4ab0edab86bd
# ╟─3a076b77-5cd6-4e10-9714-7553d2822806
# ╟─a0968f4b-b249-4488-a11b-dc109c68150f
# ╟─0e2eba31-5643-4d10-9ed4-4454ec28df12
# ╟─e0669a58-d5ac-4d01-b079-05412b413dda
# ╟─ced08ac1-8db5-42ae-b8f4-71e154ab2697
# ╟─51eb4859-20b9-4cac-bde4-ef30c6fea59d
# ╟─83851755-fe6c-4751-aa8e-3226e0fd50da
# ╟─6cc59301-fcc6-4bf4-b9ac-5d1b91560e4a
# ╟─ae93afbc-051f-47b2-a6cc-4b9f736ef9a8
# ╟─7a00bb54-553f-47f5-b5db-b40d226f4183
# ╟─fdb39284-201b-432f-bff6-986ddbc49a7d
# ╟─49faa7ea-0f22-45ca-9ab5-338d0db25564
# ╟─14db2d66-eea6-43b1-9caf-2039709d1ddb
# ╟─a2287fe8-5aa2-4259-bf7c-f715cc866243
# ╟─3c856d47-c6c2-40d3-b547-843f9654f48d
# ╟─0740f2e6-bce0-4590-acf1-ad4d7cb7c523
# ╟─ef0f33bd-8e7c-4b88-9566-df8b101828c0
# ╟─7a5a4bfa-6640-494e-88e6-1e6d5c75bb77
# ╟─5211098f-a4b9-4d1e-aaf2-d208b51a2d85
# ╟─95e1ca30-9442-4f39-9af0-34bd202fcc24
# ╟─115fa61e-8e82-42b2-8eea-9c7e21d97ea8
# ╟─f7f06be1-c8fa-4eee-953f-0d5ea26fafbf
# ╟─331846bd-406f-44e5-94db-78cf0a2594f6
# ╟─5f5e16ec-2730-4a17-bd64-a751426a033f
# ╟─e3277bb4-301a-4a1e-a838-311832b6d6aa
# ╟─85954bdb-d649-4772-a1cd-0bda5d9917e9
# ╟─8c831fdb-0bfa-4f36-b720-e82fcf5d2427
# ╟─50f1bd7f-a479-453d-a8ea-57c37a4e330c
# ╟─c26c1ea7-575f-495d-86bc-987aca991664
# ╟─13af6325-aa2c-44e2-8e2b-d05d8c5ae4a7
# ╟─3a54cad5-141f-4258-b767-49b36b6ceccc
