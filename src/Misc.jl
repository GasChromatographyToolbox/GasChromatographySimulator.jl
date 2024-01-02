# functions from GasChromatographyTools.jl:
"""
	common(s_1, s_1)

Compare two arrays and return the common elements.
"""
function common(s_1, s_2)
	common_s = String[]
	if length(s_1) >= length(s_2)
		for i=1:length(s_1)
			if s_1[i] in s_2
				push!(common_s, s_1[i])
			end
		end
	else
		for i=1:length(s_2)
			if s_2[i] in s_1
				push!(common_s, s_2[i])
			end
		end
	end
	return common_s
end

"""
	compare_peaklist(pl_1, pl_2)

Compare two peaklists (results of GasChromatographySimulator.jl) and calculate
absolute and relative differences of retention times and peak widths.
"""
function compare_peaklist(pl_1, pl_2)
	name = pl_1.Name
	tR1 = pl_1.tR
	τR1 = pl_1.τR
	tR2 = Array{Real}(undef, size(pl_1)[1])
	τR2 = Array{Real}(undef, size(pl_1)[1])
	for i=1:size(pl_1)[1]
		i2 = findfirst(name[i].==pl_2.Name)
		tR2[i] = pl_2.tR[i2]
		τR2[i] = pl_2.τR[i2]
	end
	ΔtR = tR1 .- tR2
	ΔτR = τR1 .- τR2
	rel_tR = ΔtR.*100.0./tR1
	rel_τR = ΔτR.*100.0./τR1
	compare_pl = DataFrame(Name=name, tR1=tR1, tR2=tR2, ΔtR=ΔtR, rel_tR=rel_tR, τR1=τR1, τR2=τR2, ΔτR=ΔτR, rel_τR=rel_τR)
	return compare_pl
end

"""
	compare_measurement_simulation(meas, peaklist)

Compare the retention times of measured and simulated substances.

# Arguments
* `meas`: DataFrame consisting at least of the columns `:Name` and `:RT`
  (measured retention time in s)
* `peaklist`: DataFrame as result from GasChromatographySimulator.jl with the
  columns `:Name` and `:tR` (simulated retention time in s).
  
The comparison is done by searching the same `Name` of the substance in both
DataFrames and calculating the absolute difference of the retention times (in s)
and the relative difference (in %).
"""
function compare_measurement_simulation(meas, peaklist)
	name = meas.Name
	tRm = meas.RT
	tRs = Array{Real}(undef, size(meas)[1])
	for i=1:size(meas)[1]
		i2 = findfirst(name[i].==peaklist.Name)
		if typeof(i2) == Nothing
			tRs[i] = NaN
		else
			tRs[i] = peaklist.tR[i2]
		end
	end
	compare_df = DataFrame(Name=name, measured_tR=tRm, simulated_tR=tRs, ΔtR=tRm.-tRs, rel_tR=(tRm.-tRs)./tRm.*100.0)
	return compare_df
end

"""
	load_custom_CI_database(custom_database_url)

Load a custom database for ChemicalIdentifiers.jl from the location `custom_database_url`, if the custom database is not already loaded.	
"""
function load_custom_CI_database(custom_database_url)
	#if !(:custom in keys(ChemicalIdentifiers.DATA_DB))
		ChemicalIdentifiers.load_data!(:custom, url = custom_database_url)
		ChemicalIdentifiers.load_db!(:custom)
	#end
end


function pares_No_from_sub_values(sub_values; separator=" - ")
	id = Array{Int}(undef, length(sub_values))
	for i=1:length(sub_values)
		id[i] = parse(Int, split(sub_values[i], separator)[1])
	end
	id
end