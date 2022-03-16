# functions from GasChromatographyTools.jl:
##---begin-UI-functions-----------------------------------------------------------------------------
## functions defining PlutoUI widgets for Pluto notebooks
"""
    UI_Column(sp)

Construct a combined PlutoUI widget for the settings of the GC Column with then selectable stationary phases `sp`. 
	
# UI fields
* ``L``: column length in m.
* ``d``: column diameter in mm.
* ``d_f``: film thickness in µm.
* stat. phase: stationary phase of the column
* Gas: mobile phase
"""
function UI_Column(sp; default=(10.0, 0.25, 0.25, 1, "He"))
		PlutoUI.combine() do Child
			@htl("""
			<h3>Column settings</h3>
			L [m]: $(
				Child(NumberField(0.1:0.1:100.0; default=default[1]))
			)  d [mm]: $(
				Child(NumberField(0.01:0.01:1.00; default=default[2]))
			)  d_f [µm]: $(
				Child(NumberField(0.01:0.01:1.00; default=default[3]))
			)  stat. phase: $(
				Child(Select(sp; default=sp[default[4]]))
			)  Gas: $(
				Child(Select(["He", "H2", "N2"]; default=default[5]))
			) 
			
			""")
	end
end

"""
    UI_Program(; default=("0 60 300 300 120", "40 40 170 300 300", "0 0 40 60 0", "-3 -3 -3 -3 -3", "18 18 58 98 98", "0 0 0 0 0"))

Construct a combined PlutoUI widget for the settings of the program of a GC
column with or without a thermal gradient (depending on the default tuple).

# With thermal gradient
For default as a **tuple of six strings** the folwing fields will be shown in the
widget:
* example default tupel: default=("0 60 300 300 120", "40 40 170 300 300", "0 0 40 60 0", "-3 -3 -3 -3 -3", "18 18 58 98 98", "0 0 0 0 0")
*`time steps`: the time steps after which duration the values of temperature,
inlet pressure, ΔT and α are achieved by linear interpolation (in s).
*`temperature steps`: the temperature steps (in °C). 
*`ΔT steps`: the steps of the temperature difference (in °C) between column inlet
and outlet.
*`α steps`: the steps of the gradient profile (α = 0 ... linear change of
temperature along column, α < 0 ... concave exponential profile, α > 0 ...
convexe exponential profile).
*``p_{in}`` steps: the steps of the inlet pressure (in kPa(g))
*``p_{out}`` steps: the steps of the outlet pressure (in kPa(a))

# Without thermal gradient
For a default as a **tuple of four strings** the folwing fields will be shown in the
widget:
* example default tupel: default=("0 60 600 120", "40 40 300 300", "18 18 98 98", "vacuum")
*`time steps`: the time steps after which duration the values of temperature,
inlet pressure, ΔT and α are achieved by linear interpolation (in s).
*`temperature steps`: the temperature steps (in °C). 
``p_{in}`` steps: the steps of the inlet pressure (in kPa(g))
`column outlet` selection of the outlet of the colum, "vacuum" (``p_{out} =
0.0`` kPa(a)) or "atmosphere" (``p_{out} = 101.3`` kPa(a)).

"""
function UI_Program(; default=("0 60 300 300 120", "40 40 170 300 300", "0 0 40 60 0", "-3 -3 -3 -3 -3", "18 18 58 98 98", "0 0 0 0 0"))
	if length(default) == 6
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
			) p_in steps [kPa(g)]

			$(
				Child(TextField((50,1); default=default[6]))
			) p_out steps [kPa(a)]
			</ul>	
			""")
		end
	elseif length(default) == 4
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
			) p_in steps [kPa(g)]
	
			$(
				Child(Select(["vacuum", "atmosphere"]; default=default[4]))
				) column outlet
			</ul>
			""")
		end
	end
end

"""
    UI_Substance(sol; default=(1:4, ))

Construct a combined PlutoUI widget for the settings of the substances separated
in the simulated GC system with the selectable substances `subs`. 

Depending on the tupel of `default` the widget is setup. 

For `default = (1:4, 0.0, 0.0)` the UI fields are:
* Select Substances: Selection of the substances, which will be simulated,
  default selection = 1st to 4th substance.
* Injection time: Start time (in s) of the simulation. The same for all selected
  substances. Default is 0.0 s.
* Injection width: Peak width (in s) of all selected substances at the time of
  injection. Default is 0.0 s. 

For `default = (1:4,)` the UI fields are:
* Select Substances: Selection of the substances, which will be simulated,
  default selection = 1st to 4th substance. 
"""
function UI_Substance(sol; default=(1:4,))
	if length(sol)>10
		select_size = 10
	else
		select_size = length(sol)
	end
	if length(default) == 3
		PlutoUI.combine() do Child
			@htl("""
			<h3>Substance settings</h3> 
			
			Select Substances: $(
				Child(MultiSelect(sol; default=sol[default[1]], size=select_size))
			) 
			Injection time [s]: $(
				Child(NumberField(0.0:0.1:100.0; default=default[2]))
			) and Injection width [s]: $(
				Child(NumberField(0.00:0.01:10.0; default=default[3]))
			) 
			""")
		end
	elseif length(default) == 1
		PlutoUI.combine() do Child
			@htl("""
			<h3>Substance settings</h3> 
			
			Select Substances: $(
				Child(MultiSelect(sol; default=sol[default[1]], size=select_size))
			) 
			""")
		end
	end
end

"""
    UI_Options()

Construct a combined PlutoUI widget for the settings of the options for the simulation.    
"""
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
		""")
	end
end

"""
	setting_prog(prog_values)

Translates the Program parameters from a tuple defined by a PlutoUI widget into
the structure GasChromatographySimulator.Program.
"""
function setting_prog(prog_values, L)
	# make different variants based on the size/composition of `prog_values``
	time_steps = parse.(Float64, split(prog_values[1]))
	temp_steps = parse.(Float64, split(prog_values[2]))
	pin_steps = parse.(Float64, split(prog_values[3])).*1000.0.+101300.0
	if prog_values[4] == "vacuum"
		pout_steps = zeros(length(time_steps))	
	elseif prog_values[4] == "atmosphere"
		pout_steps = 101300.0.*ones(length(time_steps))
	end
	prog = GasChromatographySimulator.Program( 	time_steps,
												temp_steps,
												pin_steps,
												pout_steps,
												L
												)
	return prog
end
##---end-UI-functions-------------------------------------------------------------------------------