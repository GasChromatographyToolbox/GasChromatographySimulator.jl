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

# ╔═╡ c9246396-3c01-4a36-bc9c-4ed72fd9e325
md"""
# Gas Chromatography Simulator

An example Simulation of a Gas Chromatography (GC) System with a thermal gradient.
"""

# ╔═╡ 8b3011fd-f3df-4ab0-b611-b943d5f3d470
md"""
## Settings
"""

# ╔═╡ 64a78d65-1553-4707-bade-10694d42daa4
sp = ["Rxi17ms", "SPB50"]

# ╔═╡ de578cc9-0662-4346-8e0b-abf2c435cfdb
function sys_set(sp)
	PlutoUI.combine() do Child
		md"""
		### System settings 
		
		``L`` [m]: $(
			Child(NumberField(0.1:0.1:100.0; default=4.0))
		) ``d`` [mm]: $(
			Child(NumberField(0.01:0.01:1.00; default=0.10))
		) ``d_f`` [µm]: $(
			Child(NumberField(0.01:0.01:1.00; default=0.10))
		) stat. phase: $(
			Child(Select(sp))
		) 

		Gas: $(
			Child(Select(["He", "H2", "N2"]))
		) 
			
		"""
	end
end

# ╔═╡ e0669a58-d5ac-4d01-b079-05412b413dda
@bind sys_values confirm(sys_set(sp))

# ╔═╡ d9549f4c-f910-452c-9a29-e8f4ba4156d1
sys = GasChromatographySimulator.System(sys_values[1], sys_values[2]*1e-3, sys_values[3]*1e-6, sys_values[4], sys_values[5])

# ╔═╡ 1ec6889b-eb37-444e-9821-0741bb0a31ca
function prog_set()
	PlutoUI.combine() do Child
		md"""
		### Program settings 
		
		$(
			Child(TextField((30,1); default="0 10 60 20"))
		) time steps [s] 
		
		$(
			Child(TextField((30,1); default="40 40 300 300"))
		) temperature steps [°C]
		
		$(
			Child(TextField((30,1); default="0 0 60 60"))
		) ``ΔT`` steps [°C]
		
		$(
			Child(TextField((30,1); default="0 0 -3 -3"))
		) ``α`` steps

		$(
			Child(TextField((30,1); default="200 200 200 200"))
		) ``p_{in}`` steps [kPa(g)]

		$(
			Child(TextField((30,1); default="101.3 101.3 101.3 101.3"))
		)``p_{out}`` steps [kPa(a)]
			
		"""
	end
end

# ╔═╡ a7e1f0ee-714e-4b97-8741-d4ab5321d5e0
@bind prog_values confirm(prog_set())

# ╔═╡ ef212663-caa3-4413-bf1b-bc4ed90e9580
prog_values
# parse the strings into Float64-Arrays
# call GasChromatographySimulator.Program

# ╔═╡ d0bf37da-dc99-4725-89f7-ac48c48b6e4f
prog_values[1]

# ╔═╡ d1791020-94af-401d-8213-c47a74d4edf6
function sub_set(sol)
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
		) Injection width [s]: $(
			Child(NumberField(0.00:0.01:10.0; default=0.0))
		) 
		"""
	end
end

# ╔═╡ d4ee5d85-a9f9-4038-a7da-f5e102dcc751
sol = ["C6", "C7", "C8", "C9", "C10", "C11", "C12"]
# load the 'Database_test.csv'
# take `sp` from the database
# take `sol` from the database (all solutes for the selected stationary phase)

# ╔═╡ 7a00bb54-553f-47f5-b5db-b40d226f4183
@bind sub_values confirm(sub_set(sol))

# ╔═╡ 84171b9f-1c0f-4250-a8cd-67fd3f646ed9
sub_values
# call GasChromatographySimulator.load_solute_database()

# ╔═╡ d33b3294-1293-47c1-aa39-f7dc6025986d
function opt_set()
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

# ╔═╡ 3e053ac1-db7b-47c1-b52c-00e26b59912f
@bind opt_values confirm(opt_set())

# ╔═╡ 69a1522a-c6aa-4a0d-baf6-4b9dec0a466a
opt = GasChromatographySimulator.Options(;abstol=10.0^opt_values[1], reltol=10.0^opt_values[2], Tcontrol=opt_values[3])

# ╔═╡ Cell order:
# ╠═7272450e-73b1-11ec-080d-1d1efd32e836
# ╟─c9246396-3c01-4a36-bc9c-4ed72fd9e325
# ╟─8b3011fd-f3df-4ab0-b611-b943d5f3d470
# ╠═64a78d65-1553-4707-bade-10694d42daa4
# ╠═de578cc9-0662-4346-8e0b-abf2c435cfdb
# ╠═e0669a58-d5ac-4d01-b079-05412b413dda
# ╠═d9549f4c-f910-452c-9a29-e8f4ba4156d1
# ╠═1ec6889b-eb37-444e-9821-0741bb0a31ca
# ╠═a7e1f0ee-714e-4b97-8741-d4ab5321d5e0
# ╠═ef212663-caa3-4413-bf1b-bc4ed90e9580
# ╠═d0bf37da-dc99-4725-89f7-ac48c48b6e4f
# ╠═d1791020-94af-401d-8213-c47a74d4edf6
# ╠═d4ee5d85-a9f9-4038-a7da-f5e102dcc751
# ╠═7a00bb54-553f-47f5-b5db-b40d226f4183
# ╠═84171b9f-1c0f-4250-a8cd-67fd3f646ed9
# ╠═d33b3294-1293-47c1-aa39-f7dc6025986d
# ╠═3e053ac1-db7b-47c1-b52c-00e26b59912f
# ╠═69a1522a-c6aa-4a0d-baf6-4b9dec0a466a
