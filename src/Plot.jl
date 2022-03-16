##---begin-Plot-functions---------------------------------------------------------------------------
"""
    chromatogram(t, tR, τR)

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
    plot_chromatogram(peaklist, tlims; annotation=true, number=true, mirror=false, offset=0.0)

Plot the chromatogram of the peaks listed in `peaklist` over the time tupel `tlims = (t_start, t_end`). 

# Arguments
* `peaklist`: DataFrame with the names, retention times and peak widths of the simulated substances.
* `tlims`: Tuple defining the start and end time of the plotted chromatogram.
* `annotation`: Boolean, switching the annotation of the peaks on/off; default = true.
* `number`: Boolean, switching the type of the annotation between _number_ of the substance in the peaklist (`number = true`) or _name_ of the substance (`number = false`); default = true.
* `mirror`: Boolean, if `mirror = true` the chromatogram is multiplied by `-1`; default = false.
* `offset`: Float64, this value is added to the chromatogram; default = 0.0.

# Output
Tupel `(p_chrom, t, chrom)`
* `p_chrom`: the plot of the chromatogram `chrom` over time `t`
* `t`: Array of time of the chromatogram
* `chrom`: Array of the abundance values of the chromatogram
"""
function plot_chromatogram(peaklist, tlims; annotation=true, number=true, mirror=false, offset=0.0)
	t₀ = tlims[1]
	tMax = tlims[2]
	t = t₀:tMax/10000:tMax
	if mirror==true
		sgn = -1.0
	else
		sgn = 1.0
	end
	chrom = sgn.*chromatogram(collect(t), peaklist.tR, peaklist.τR) .+ offset
	chrom_tR = sgn.*chromatogram(peaklist.tR, peaklist.tR, peaklist.τR) .+ offset
	p_chrom = plot(t, chrom, xlabel="time in s", ylabel="abundance", legend=false)
	if annotation==true && number==true
		plot!(p_chrom, annotations = [(peaklist.tR[i], chrom_tR[i], text(i, 10, rotation=0, :center)) for i in 1:length(peaklist.tR)])
	elseif annotation==true && number==false
		plot!(p_chrom, annotations = [(peaklist.tR[i], chrom_tR[i], text(peaklist.Name[i], 10, rotation=90, :center)) for i in 1:length(peaklist.tR)])
	end
	return p_chrom, t, chrom
end

"""
    plot_chromatogram!(p_chrom, peaklist, tlims; annotation=true, number=true, mirror=false, offset=0.0)

Add the chromatogram of the peaks listed in `peaklist` over the time tupel `tlims = (t_start, t_end`) to the plot `p_chrom`. 

# Arguments
* `p_chrom`: Plot of an existing chromatogram.
* `peaklist`: DataFrame with the names, retention times and peak widths of the simulated substances.
* `tlims`: Tuple defining the start and end time of the plotted chromatogram.
* `annotation`: Boolean, switching the annotation of the peaks on/off; default = true.
* `number`: Boolean, switching the type of the annotation between _number_ of the substance in the peaklist (`number = true`) or _name_ of the substance (`number = false`); default = true.
* `mirror`: Boolean, if `mirror = true` the chromatogram is multiplied by `-1`; default = false.
* `offset`: Float64, this value is added to the chromatogram; default = 0.0.

# Output
Tupel `(t, chrom)`
* `t`: Array of time of the chromatogram
* `chrom`: Array of the abundance values of the chromatogram
"""
function plot_chromatogram!(p_chrom, peaklist, tlims; t₀=0.0, annotation=true, number=true, mirror=false, offset=0.0)
	t₀ = tlims[1]
	tMax = tlims[2]
	t = t₀:tMax/10000:tMax
	if mirror==true
		sgn = -1.0
	else
		sgn = 1.0
	end
	chrom = sgn.*chromatogram(collect(t), peaklist.tR, peaklist.τR) .+ offset
	chrom_tR = sgn.*chromatogram(peaklist.tR, peaklist.tR, peaklist.τR) .+ offset
	plot!(p_chrom, t, chrom)
	if annotation==true && number==true
		plot!(p_chrom, annotations = [(peaklist.tR[i], chrom_tR[i], text(i, 10, rotation=0, :center)) for i in 1:length(peaklist.tR)])
	elseif annotation==true && number==false
		plot!(p_chrom, annotations = [(peaklist.tR[i], chrom_tR[i], text(peaklist.Name[i], 10, rotation=90, :center)) for i in 1:length(peaklist.tR)])
	end
	return t, chrom
end

"""
	plot_flow(par)

Calculate and plot the flow (in mL/min, normalized) of the carrier gas in a GC Column with a program defined in the parameters `par::GasChromatography.Parameters`.
"""
function plot_flow(par)
	t = 0.0:sum(par.prog.time_steps)/1000.0:sum(par.prog.time_steps)
	F = Array{Float64}(undef, length(t))
	for i=1:length(t)
		F[i] = flow(t[i], par.prog.T_itp, par.prog.pin_itp, par.prog.pout_itp, par.col.L, par.col.d, par.col.gas; vis=par.opt.vis)
	end
	p_flow = plot(t, F.*60e6, xlabel="time in s", ylabel="column flow in mL/min", legend=false)
	return p_flow
end

"""
	plot_pressure(prog)

Plot the inlet and outlet pressure over time of the program `prog::GasChromatographySimulator.Program`.
"""
function plot_pressure(prog)
	t = 0.0:sum(prog.time_steps)/1000.0:sum(prog.time_steps)
	pin = Array{Float64}(undef, length(t))
	pout = Array{Float64}(undef, length(t))
	for i=1:length(t)
		pin[i] = prog.pin_itp(t[i])
		pout[i] = prog.pout_itp(t[i])
	end
	p_press = plot(t, pin, xlabel="time in s", ylabel="pressure in Pa", label="inlet", legend=:right)
	plot!(p_press, t, pout, label="outlet")
	return p_press
end

"""
	plot_temperature(par; selector="T(t)")

Plot the temperature program of the GC Column. 

# Arguments
* `par::GasChromatographySimulator.Parameters`: parameters of the GC system
* `selector::String`: selection of the plot as:
	* `selector = "T(t)"`: 2D-plot of temperature `T` over time `t` at inlet (`x=0`) and outlet (`x=L`)
	* `selector = "T(x)"`: 2D-plot of temperature `T` over column position `x` at the `time_steps` of the program
	* `selector = "T(x,t)"`: 3D-plot of temperature `T` over column position `x` and `t`
"""
function plot_temperature(par::GasChromatographySimulator.Parameters; selector="T(t)")
	if selector=="T(x)"
		nx = 0.0:par.col.L/1000:par.col.L
		p_temp = plot(xlabel="x in m", ylabel="T in °C", legend=:top)
		for i=1:length(par.prog.time_steps)
			T = par.prog.T_itp.(nx, cumsum(par.prog.time_steps)[i]).-273.15
			plot!(p_temp, nx, T, label="t=$(cumsum(par.prog.time_steps)[i])s")
		end
	elseif selector=="T(t)"
		nt = 0.0:sum(par.prog.time_steps)/1000:sum(par.prog.time_steps)
		T0 = par.prog.T_itp.(0.0, nt).-273.15
		p_temp = plot(nt, T0, xlabel="t in s", ylabel="T in °C", legend=:top, label="inlet", c=:red)
		TL = par.prog.T_itp.(par.col.L, nt).-273.15
		plot!(p_temp, nt, TL, label="outlet", c=:blue)
	elseif selector=="T(x,t)"
		nx = 0.0:par.col.L/1000:par.col.L
		nt = 0.0:sum(par.prog.time_steps)/1000:sum(par.prog.time_steps)
		T = Array{Float64}(undef, length(nx), length(nt))
		for j=1:length(nt)
			for i=1:length(nx)
				T[i,j] = par.prog.T_itp(nx[i], nt[j])-273.15
			end
		end
		p_temp = plot(nx, nt, T', st=:surface, xlabel="x in m", ylabel="t in s", zlabel="T in °C")
	end
	return p_temp
end

# functions from GasChromatographyTools.jl:
"""
	local_plots(xx, yy, sol, par)

Show additional 'local' plots of selected `yy` quantities over selected `xx`
quantities.

# Arguments
* `xx`: Selected quantity shown on the x-axis. Possible values: "z", "t", "T",
  "τ", "σ" and "u".
* `yy`: Selected quantity shown on the y-axis. Possible values: "z", "t", "T",
  "τ", "σ" and "u".
* `sol`: The solution of the simulation.
* `par`: The parameters of the simulated GC-system.
"""   
function local_plots(xx, yy, sol, par)
	n = size(sol)[1]

	df_sol = GasChromatographySimulator.sol_extraction(sol, par)
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

"""
	velocity(df_sol, i, par)

Calculate the velocity (in m/s) coressponding to solution of the `i-th` sunstance of a
GC-system defined by `par`.
""" 
function velocity(df_sol, i, par)
	x = df_sol.z[i]
	t = df_sol.t[i]
	T_itp = par.prog.T_itp
	pin_itp = par.prog.pin_itp
	pout_itp = par.prog.pout_itp
	L = par.col.L
	d = par.col.d
	df = par.col.df
	gas = par.col.gas
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