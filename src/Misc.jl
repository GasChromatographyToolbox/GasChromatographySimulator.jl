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
	# check this function
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

"""
	isinteger(x)

Wrapper function for `isa(x, Int)`.
"""
isinteger(x) = isa(x, Int)

"""
	same_length(a, b)

Changes the length of vector `a` to the length of vector `b`. If `a` is shorter than `b` then additional zeros will be added to `a`. If `a` is longer than `b` then only the first elements `a` will be used.
"""
function same_length(a, b)
	a = if length(b) > length(a)
		[a; zeros(length(b)-length(a))]
	elseif length(b) < length(a)
		a[1:length(b)]
	else
		a
	end
end

"""
    deduplicate_knots!(x; move_knots=true)

Remove duplicate consecutive values from array x, optionally moving knots to avoid duplicates.
Similar to Interpolations.deduplicate_knots! but simpler.

# Arguments
* `x` ... array to deduplicate (modified in place)
* `move_knots` ... if true, slightly adjust duplicate values to make them unique

# Output
* Indices of unique values (returns 1:length(x) after modification)

# Examples
```julia
x = [1.0, 2.0, 2.0, 3.0]
unique_indices = deduplicate_knots!(x; move_knots=true)
# x is now [1.0, 2.0, 2.0000000001, 3.0] (approximately)
# unique_indices = 1:4
```
"""
function deduplicate_knots!(x; move_knots=true)
    if length(x) <= 1
        return 1:length(x)
    end
    
    # Modify duplicates to make them unique (similar to Interpolations.deduplicate_knots!)
    # When move_knots=true, slightly adjust duplicate values
    if move_knots
        i = 1
        while i < length(x)
            # Find the start of a duplicate sequence
            if x[i+1] == x[i]
                # Found a duplicate, find the end of the duplicate sequence
                dup_start = i
                dup_end = i + 1
                while dup_end < length(x) && x[dup_end + 1] == x[dup_start]
                    dup_end += 1
                end
                
                # Now modify all duplicates in the sequence
                if dup_end < length(x)
                    # Not at the end: adjust based on next value
                    next_val = x[dup_end + 1]
                    for j in (dup_start + 1):dup_end
                        # Use a small increment that increases with position in duplicate sequence
                        offset = 1e-10 * (next_val - x[dup_start]) * (j - dup_start) / (dup_end - dup_start + 1)
                        # If next value is same, use absolute increment
                        if next_val == x[dup_start]
                            offset = 1e-10 * abs(x[dup_start]) * (j - dup_start) + 1e-15 * (j - dup_start)
                        end
                        x[j] = x[dup_start] + offset
                    end
                else
                    # At the end: adjust based on previous difference
                    # Find the last different value
                    prev_val = x[dup_start]
                    prev_diff = 0.0
                    for j in (dup_start - 1):-1:1
                        if x[j] != x[dup_start]
                            prev_val = x[j]
                            prev_diff = x[dup_start] - x[j]
                            break
                        end
                    end
                    
                    if prev_diff == 0.0
                        # All previous values are the same, use small increment
                        for j in (dup_start + 1):dup_end
                            x[j] = x[dup_start] + 1e-10 * abs(x[dup_start]) * (j - dup_start) + 1e-15 * (j - dup_start)
                        end
                    else
                        # Use previous difference as direction
                        for j in (dup_start + 1):dup_end
                            x[j] = x[dup_start] + 1e-10 * prev_diff * (j - dup_start)
                        end
                    end
                end
                
                # Skip to after the duplicate sequence
                i = dup_end + 1
            else
                i += 1
            end
        end
    end
    
    # Return all indices (after modification, all values should be unique if move_knots=true)
    return 1:length(x)
end