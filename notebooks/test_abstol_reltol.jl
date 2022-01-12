### A Pluto.jl notebook ###
# v0.17.5

using Markdown
using InteractiveUtils

# ╔═╡ 6f9a94bd-dd28-4110-a7ca-0ed84e9c7c3f
begin
	import Pkg
    # activate the shared project environment
    Pkg.activate(Base.current_project())
	using Plots
	using GasChromatographySimulator
	using Statistics
	using Plots
	using PlutoUI
	TableOfContents()
end

# ╔═╡ 0b51a792-d9b8-4c29-928f-57aaf41f1c20
plotly()

# ╔═╡ 93ba6afc-4e9a-11ec-08d9-81c0a9bc502e
md"""
# Test of the influence of abstol and reltol on the results of the GC simulation
"""

# ╔═╡ 51ec0223-199a-4a25-8768-cd0b7e9f864d
begin
	# stretch of the GC program by factor n. Taken from GasChromatographyTools
	function stretched_program(n::Float64, par::GasChromatographySimulator.Parameters)
		# stretch the temperature program in 'par' by a factor 'n'
		if isa(par, Array)==true
			error("Select an element of the array of GC-system parameters.")
		else
			new_tsteps = n.*par.prog.time_steps
			new_T_itp = GasChromatographySimulator.temperature_interpolation(new_tsteps, par.prog.temp_steps, par.prog.gf, par.sys.L)
			new_pin_itp = GasChromatographySimulator.pressure_interpolation(new_tsteps, par.prog.pin_steps)
			new_pout_itp = GasChromatographySimulator.pressure_interpolation(new_tsteps, par.prog.pout_steps)
			new_prog = GasChromatographySimulator.Program(new_tsteps, par.prog.temp_steps, par.prog.pin_steps, par.prog.pout_steps, par.prog.gf, par.prog.a_gf, new_T_itp, new_pin_itp, new_pout_itp)
			new_par = GasChromatographySimulator.Parameters(par.sys, new_prog, par.sub, par.opt)
			return new_par
		end	
	end
end

# ╔═╡ 5e642f09-9a8f-4cca-b61f-b27c8433a2e5
begin
	opt = GasChromatographySimulator.Options(OwrenZen5(), 1e-6, 1e-4, "inlet", true)
	L = 4.0
	d = 0.1e-3
	df = 0.1e-6
	sp = "Rxi17SilMS" # ["Rxi17SilMS" -> ERR, "SLB5ms", "SPB50", "Wax", "DB5ms", "Rxi5MS", "genericLB", "genericJL"]
	gas = "He"
	sys = GasChromatographySimulator.constructor_System(L, d, df, sp, gas)

	db_path = "/Users/janleppert/Documents/GitHub/Publication_GCsim/data/Databases/" 
	db_file = "Database_append.csv"
	first_alkane = "C8"
	last_alkane = "C15"
	sub = GasChromatographySimulator.load_solute_database(db_path, db_file, sp, gas, [first_alkane, last_alkane], zeros(2), zeros(2))
	
	Tst = 273.15
	Tref = 150.0 + Tst
	pn = 101300.0
	pin = 300000.0 + pn
	pout = 0.0
	dimless_rate = 0.4
	Theat = 1000.0
	Tshift = 40.0
	tMref = GasChromatographySimulator.holdup_time(Tref, pin, pout, L, d, gas)
	rate = dimless_rate*30/tMref
	theat = Theat/rate
	Tstart = sub[1].Tchar-Tst-Tshift
	
	ΔT = 135.0
	α = 9.0
	prog0 = GasChromatographySimulator.constructor_Program([0.0, theat],[Tstart, Tstart+Theat], [pin, pin],[pout, pout],[ΔT, ΔT], [0.0, 0.0], [L, L], [α, α], opt.Tcontrol, L)
	
	par0 = GasChromatographySimulator.Parameters(sys, prog0, sub, opt)

	tR_lock = 12*tMref
	tR_tol = 1e-3
	md"""
	## Settings
	The following settings produce the problem:
	"""
end

# ╔═╡ ae803a26-c436-497a-b245-a0afec92e46f
md"""
Make simulations for a range of stretch factors ``n``.
"""

# ╔═╡ 08a0972a-55ff-42c9-838b-1a649afe9a46
md"""
## Decreased relative tolerance

The problem is not originated in the RT_lock algorithm but in the simulation itself. For certain settings a small change (like small strectch in the program) can result in bigger changes in retention time.
"""

# ╔═╡ 2d7da55e-2711-44a4-b8b9-bda6321a4c48
begin
	# decrease the relative tolerance
	opt_6_3 = GasChromatographySimulator.Options(OwrenZen5(), 1e-6, 1e-3, "inlet", true)
	par_6_3 = GasChromatographySimulator.Parameters(sys, prog0, sub, opt_6_3)
	# repeat the simulation from above
	n_6_3 = sort!(rand(0.499625:0.000000001:0.499627, 100))
	tR_6_3 = Array{Float64}(undef, length(n_6_3))
	for i=1:length(n_6_3)
	    pars = stretched_program(n_6_3[i], par_6_3)
	    sols = GasChromatographySimulator.solving_odesystem_r(pars.sys, pars.prog, pars.sub[2], pars.opt)
	    tR_6_3[i] = sols.u[end][1]
	end
end

# ╔═╡ 7f968c24-1ddb-42f7-83ef-403629615ea3
p1 = plot(n_6_3, tR_6_3, markers=:o, label="6_3")

# ╔═╡ 1639930a-6f3f-43a1-a750-b36b66f9541f
begin
	# decrease the relative tolerance
	opt_6_4 = GasChromatographySimulator.Options(OwrenZen5(), 1e-6, 1e-4, "inlet", true)
	par_6_4 = GasChromatographySimulator.Parameters(sys, prog0, sub, opt_6_4)
	# repeat the simulation from above
	n_6_4 = sort!(rand(0.499625:0.000000001:0.499627, 100))
	tR_6_4 = Array{Float64}(undef, length(n_6_4))
	for i=1:length(n_6_4)
	    pars = stretched_program(n_6_4[i], par_6_4)
	    sols = GasChromatographySimulator.solving_odesystem_r(pars.sys, pars.prog, pars.sub[2], pars.opt)
	    tR_6_4[i] = sols.u[end][1]
	end
end

# ╔═╡ 1de2a221-d047-4206-9a50-d98274e995cf
plot!(p1, n_6_4, tR_6_4, markers=:o, label="6_4")

# ╔═╡ 38349c2f-27eb-439c-b4fd-3306e64994e7
begin
	# decrease the relative tolerance
	opt_6_5 = GasChromatographySimulator.Options(OwrenZen5(), 1e-6, 1e-5, "inlet", true)
	par_6_5 = GasChromatographySimulator.Parameters(sys, prog0, sub, opt_6_5)
	# repeat the simulation from above
	n_6_5 = sort!(rand(0.499625:0.000000001:0.499627, 100))
	tR_6_5 = Array{Float64}(undef, length(n_6_5))
	for i=1:length(n_6_5)
	    pars = stretched_program(n_6_5[i], par_6_5)
	    sols = GasChromatographySimulator.solving_odesystem_r(pars.sys, pars.prog, pars.sub[2], pars.opt)
	    tR_6_5[i] = sols.u[end][1]
	end
end

# ╔═╡ 67097337-b99d-4263-b39e-fe142dc71ccd
plot!(p1, n_6_5, tR_6_5, markers=:o, label="6_5")

# ╔═╡ 41b8302c-0859-4aa2-a149-dc1037b0aaa8
begin
	# decrease the relative tolerance
	opt_6_6 = GasChromatographySimulator.Options(OwrenZen5(), 1e-6, 1e-6, "inlet", true)
	par_6_6 = GasChromatographySimulator.Parameters(sys, prog0, sub, opt_6_6)
	# repeat the simulation from above
	n_6_6 = sort!(rand(0.499625:0.000000001:0.499627, 100))
	tR_6_6 = Array{Float64}(undef, length(n_6_6))
	for i=1:length(n_6_6)
	    pars = stretched_program(n_6_6[i], par_6_6)
	    sols = GasChromatographySimulator.solving_odesystem_r(pars.sys, pars.prog, pars.sub[2], pars.opt)
	    tR_6_6[i] = sols.u[end][1]
	end
end

# ╔═╡ c4d35d83-792d-4bb1-a3a8-6771c194fd6f
plot!(p1, n_6_6, tR_6_6, markers=:o, label="6_6")

# ╔═╡ 224ddaae-480e-456d-a120-132dbf6b0966
begin
	# decrease the relative tolerance
	opt_7_4 = GasChromatographySimulator.Options(OwrenZen5(), 1e-7, 1e-4, "inlet", true)
	par_7_4 = GasChromatographySimulator.Parameters(sys, prog0, sub, opt_7_4)
	# repeat the simulation from above
	n_7_4 = sort!(rand(0.499625:0.000000001:0.499627, 100))
	tR_7_4 = Array{Float64}(undef, length(n_7_4))
	for i=1:length(n_7_4)
	    pars = stretched_program(n_7_4[i], par_7_4)
	    sols = GasChromatographySimulator.solving_odesystem_r(pars.sys, pars.prog, pars.sub[2], pars.opt)
	    tR_7_4[i] = sols.u[end][1]
	end
end

# ╔═╡ b695f569-2701-4307-a537-791f1c17e515
begin
	p2 = plot(n_6_4, tR_6_4, markers=:o, label="6_4")
	plot!(p2, n_7_4, tR_7_4, markers=:o, label="7_4")
	p2
end

# ╔═╡ 78f6adb7-8f7b-4e9d-ab27-70b3975e6cc3
begin
	# decrease the relative tolerance
	opt_8_4 = GasChromatographySimulator.Options(OwrenZen5(), 1e-8, 1e-4, "inlet", true)
	par_8_4 = GasChromatographySimulator.Parameters(sys, prog0, sub, opt_8_4)
	# repeat the simulation from above
	n_8_4 = sort!(rand(0.499625:0.000000001:0.499627, 100))
	tR_8_4 = Array{Float64}(undef, length(n_8_4))
	for i=1:length(n_8_4)
	    pars = stretched_program(n_8_4[i], par_8_4)
	    sols = GasChromatographySimulator.solving_odesystem_r(pars.sys, pars.prog, pars.sub[2], pars.opt)
	    tR_8_4[i] = sols.u[end][1]
	end
end

# ╔═╡ 275a306b-4176-4392-a08f-0be88134a908
plot!(p2, n_8_4, tR_8_4, markers=:o, label="8_4")

# ╔═╡ 03eb99d3-6fcc-4e33-b84f-e9ce4593c0f5
begin
	# decrease the relative tolerance
	opt_9_4 = GasChromatographySimulator.Options(OwrenZen5(), 1e-9, 1e-4, "inlet", true)
	par_9_4 = GasChromatographySimulator.Parameters(sys, prog0, sub, opt_9_4)
	# repeat the simulation from above
	n_9_4 = sort!(rand(0.499625:0.000000001:0.499627, 100))
	tR_9_4 = Array{Float64}(undef, length(n_9_4))
	for i=1:length(n_9_4)
	    pars = stretched_program(n_9_4[i], par_9_4)
	    sols = GasChromatographySimulator.solving_odesystem_r(pars.sys, pars.prog, pars.sub[2], pars.opt)
	    tR_9_4[i] = sols.u[end][1]
	end
end

# ╔═╡ c989e68b-d661-49d2-bfc1-2872b9b09242
plot!(p2, n_9_4, tR_9_4, markers=:o, label="9_4")

# ╔═╡ acbd4b5a-9ae3-460c-9dac-412675a3721f
begin
	# decrease the relative tolerance
	opt_10_4 = GasChromatographySimulator.Options(OwrenZen5(), 1e-10, 1e-4, "inlet", true)
	par_10_4 = GasChromatographySimulator.Parameters(sys, prog0, sub, opt_10_4)
	# repeat the simulation from above
	n_10_4 = sort!(rand(0.499625:0.000000001:0.499627, 100))
	tR_10_4 = Array{Float64}(undef, length(n_10_4))
	for i=1:length(n_10_4)
	    pars = stretched_program(n_10_4[i], par_10_4)
	    sols = GasChromatographySimulator.solving_odesystem_r(pars.sys, pars.prog, pars.sub[2], pars.opt)
	    tR_10_4[i] = sols.u[end][1]
	end
end

# ╔═╡ b8497395-6926-4c8f-94e9-d51bbcaf8cfd
plot!(p2, n_10_4, tR_10_4, markers=:o, label="10_4")

# ╔═╡ 6cf71fb9-57ce-45f8-8e1f-61ef51788bf1
begin
	# decrease the relative tolerance
	opt_9_3 = GasChromatographySimulator.Options(OwrenZen5(), 1e-9, 1e-3, "inlet", true)
	par_9_3 = GasChromatographySimulator.Parameters(sys, prog0, sub, opt_9_3)
	# repeat the simulation from above
	n_9_3 = sort!(rand(0.499625:0.000000001:0.499627, 100))
	tR_9_3 = Array{Float64}(undef, length(n_9_3))
	for i=1:length(n_9_3)
	    pars = stretched_program(n_9_3[i], par_9_3)
	    sols = GasChromatographySimulator.solving_odesystem_r(pars.sys, pars.prog, pars.sub[2], pars.opt)
	    tR_9_3[i] = sols.u[end][1]
	end
end

# ╔═╡ 60819ce0-5d84-4def-851f-3e9abad32204
begin
	p3 = plot(n_9_3, tR_9_3, markers=:o, label="9_3")
	plot!(p3, n_9_4, tR_9_4, markers=:o, label="9_4")
	p3
end

# ╔═╡ f7fc6573-97e4-4c0a-8713-d436c43a086d
begin
	# decrease the relative tolerance
	opt_9_5 = GasChromatographySimulator.Options(OwrenZen5(), 1e-9, 1e-5, "inlet", true)
	par_9_5 = GasChromatographySimulator.Parameters(sys, prog0, sub, opt_9_5)
	# repeat the simulation from above
	n_9_5 = sort!(rand(0.499625:0.000000001:0.499627, 100))
	tR_9_5 = Array{Float64}(undef, length(n_9_5))
	for i=1:length(n_9_5)
	    pars = stretched_program(n_9_5[i], par_9_5)
	    sols = GasChromatographySimulator.solving_odesystem_r(pars.sys, pars.prog, pars.sub[2], pars.opt)
	    tR_9_5[i] = sols.u[end][1]
	end
end

# ╔═╡ bc9c3068-37cb-44a8-8240-98584ca0a1c6
plot!(p3, n_9_5, tR_9_5, markers=:o, label="9_5")

# ╔═╡ d9466a04-fa39-4d4b-a98f-b7ab2b0eca1a
begin
	# decrease the relative tolerance
	opt_9_6 = GasChromatographySimulator.Options(OwrenZen5(), 1e-9, 1e-6, "inlet", true)
	par_9_6 = GasChromatographySimulator.Parameters(sys, prog0, sub, opt_9_6)
	# repeat the simulation from above
	n_9_6 = sort!(rand(0.499625:0.000000001:0.499627, 100))
	tR_9_6 = Array{Float64}(undef, length(n_9_6))
	for i=1:length(n_9_6)
	    pars = stretched_program(n_9_6[i], par_9_6)
	    sols = GasChromatographySimulator.solving_odesystem_r(pars.sys, pars.prog, pars.sub[2], pars.opt)
	    tR_9_6[i] = sols.u[end][1]
	end
end

# ╔═╡ 3faa131d-c6f1-4e2f-ab23-a535a1f80a99
plot!(p3, n_9_6, tR_9_6, markers=:o, label="9_6")

# ╔═╡ c4dd70b2-95c0-4a1b-aa31-ea5a3e63271a
begin
	# decrease the relative tolerance
	opt_9_7 = GasChromatographySimulator.Options(OwrenZen5(), 1e-9, 1e-7, "inlet", true)
	par_9_7 = GasChromatographySimulator.Parameters(sys, prog0, sub, opt_9_7)
	# repeat the simulation from above
	n_9_7 = sort!(rand(0.499625:0.000000001:0.499627, 100))
	tR_9_7 = Array{Float64}(undef, length(n_9_7))
	for i=1:length(n_9_7)
	    pars = stretched_program(n_9_7[i], par_9_7)
	    sols = GasChromatographySimulator.solving_odesystem_r(pars.sys, pars.prog, pars.sub[2], pars.opt)
	    tR_9_7[i] = sols.u[end][1]
	end
end

# ╔═╡ bc0bc452-e409-472e-a7cd-dc61256c6121
plot!(p3, n_9_7, tR_9_7, markers=:o, label="9_7")

# ╔═╡ 0aeb3ece-93a4-4638-bd90-1fdf13231e2c
md"""
**This is a good test inviroment for variations of abstol and reltol.**

It seems, that, if the result is 'noisy' (which is not always the case), than we have two branches of results for retention times for variing program stretch factors n. 

-> new notebook in GasChromatographySimulator.jl
"""

# ╔═╡ a0e85c4b-7332-42a6-8fa0-77a9bb0c89a8
md"""
# Conclusion
"""

# ╔═╡ 42a1fccd-e782-42fb-91de-18c8b33d507d
md"""
# End
"""

# ╔═╡ Cell order:
# ╠═6f9a94bd-dd28-4110-a7ca-0ed84e9c7c3f
# ╠═0b51a792-d9b8-4c29-928f-57aaf41f1c20
# ╠═93ba6afc-4e9a-11ec-08d9-81c0a9bc502e
# ╠═51ec0223-199a-4a25-8768-cd0b7e9f864d
# ╟─5e642f09-9a8f-4cca-b61f-b27c8433a2e5
# ╠═ae803a26-c436-497a-b245-a0afec92e46f
# ╟─08a0972a-55ff-42c9-838b-1a649afe9a46
# ╠═2d7da55e-2711-44a4-b8b9-bda6321a4c48
# ╠═7f968c24-1ddb-42f7-83ef-403629615ea3
# ╠═1639930a-6f3f-43a1-a750-b36b66f9541f
# ╠═1de2a221-d047-4206-9a50-d98274e995cf
# ╠═38349c2f-27eb-439c-b4fd-3306e64994e7
# ╠═67097337-b99d-4263-b39e-fe142dc71ccd
# ╠═41b8302c-0859-4aa2-a149-dc1037b0aaa8
# ╠═c4d35d83-792d-4bb1-a3a8-6771c194fd6f
# ╠═224ddaae-480e-456d-a120-132dbf6b0966
# ╠═b695f569-2701-4307-a537-791f1c17e515
# ╠═78f6adb7-8f7b-4e9d-ab27-70b3975e6cc3
# ╠═275a306b-4176-4392-a08f-0be88134a908
# ╠═03eb99d3-6fcc-4e33-b84f-e9ce4593c0f5
# ╠═c989e68b-d661-49d2-bfc1-2872b9b09242
# ╠═acbd4b5a-9ae3-460c-9dac-412675a3721f
# ╠═b8497395-6926-4c8f-94e9-d51bbcaf8cfd
# ╠═6cf71fb9-57ce-45f8-8e1f-61ef51788bf1
# ╠═60819ce0-5d84-4def-851f-3e9abad32204
# ╠═f7fc6573-97e4-4c0a-8713-d436c43a086d
# ╠═bc9c3068-37cb-44a8-8240-98584ca0a1c6
# ╠═d9466a04-fa39-4d4b-a98f-b7ab2b0eca1a
# ╠═3faa131d-c6f1-4e2f-ab23-a535a1f80a99
# ╠═c4dd70b2-95c0-4a1b-aa31-ea5a3e63271a
# ╠═bc0bc452-e409-472e-a7cd-dc61256c6121
# ╠═0aeb3ece-93a4-4638-bd90-1fdf13231e2c
# ╠═a0e85c4b-7332-42a6-8fa0-77a9bb0c89a8
# ╠═42a1fccd-e782-42fb-91de-18c8b33d507d
