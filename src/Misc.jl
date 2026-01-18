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

#---Begin-Simple-Linear-Interpolation-Implementation---
"""
    SimpleLinearInterpolation

A simple linear interpolation struct that is fully compatible with ForwardDiff.jl.
Supports 1D and 2D interpolation with flat extrapolation.

# Fields
* `grid`: Tuple of grid points for each dimension (Vector for 1D, Tuple of 2 Vectors for 2D)
* `values`: Array of values at grid points (Vector for 1D, Matrix for 2D)
"""
struct SimpleLinearInterpolation{D, T}
    grid::NTuple{D, Vector{T}}
    values::Array{T, D}
end

# Equality comparison for SimpleLinearInterpolation
function Base.:(==)(itp1::SimpleLinearInterpolation, itp2::SimpleLinearInterpolation)
    # Check dimensions match
    if typeof(itp1).parameters[1] != typeof(itp2).parameters[1]
        return false
    end
    D = typeof(itp1).parameters[1]
    
    # Compare grids
    for i in 1:D
        if length(itp1.grid[i]) != length(itp2.grid[i])
            return false
        end
        if !isapprox(itp1.grid[i], itp2.grid[i], rtol=1e-14, atol=1e-14)
            return false
        end
    end
    
    # Compare values
    if size(itp1.values) != size(itp2.values)
        return false
    end
    return isapprox(itp1.values, itp2.values, rtol=1e-14, atol=1e-14)
end

# Helper function to extract numeric value for comparisons (works with Dual numbers and Measurements)
# Order matters: more specific types first
_extract_value(x::Measurements.Measurement) = Measurements.value(x)
_extract_value(x::Real) = x
_extract_value(x) = ForwardDiff.value(x)

# 1D linear interpolation
function (itp::SimpleLinearInterpolation{1, T})(t) where T
    grid = itp.grid[1]
    values = itp.values
    
    # Extract numeric value for index finding
    t_val = _extract_value(t)
    
    # Flat extrapolation: return boundary values when out of bounds
    if t_val <= grid[1]
        return values[1]
    elseif t_val >= grid[end]
        return values[end]
    end
    
    # Find the interval containing t
    # Use searchsortedlast to find the largest index where grid[i] <= t_val
    i = searchsortedlast(grid, t_val)
    
    # Handle edge case where t_val equals the last grid point
    if i == length(grid)
        return values[end]
    end
    
    # Linear interpolation (preserves Dual number structure)
    t1, t2 = grid[i], grid[i+1]
    v1, v2 = values[i], values[i+1]
    
    # Compute weight: w = (t - t1) / (t2 - t1)
    # This arithmetic preserves ForwardDiff.Dual structure
    w = (t - t1) / (t2 - t1)
    
    # Interpolate: v = v1 + w * (v2 - v1)
    return v1 + w * (v2 - v1)
end

# 2D bilinear interpolation
function (itp::SimpleLinearInterpolation{2, T})(x, t) where T
    x_grid, t_grid = itp.grid
    values = itp.values
    
    # Extract numeric values for index finding
    x_val = _extract_value(x)
    t_val = _extract_value(t)
    
    # Flat extrapolation: clamp to boundaries
    x_val_clamped = clamp(x_val, x_grid[1], x_grid[end])
    t_val_clamped = clamp(t_val, t_grid[1], t_grid[end])
    
    # Find intervals
    i = searchsortedlast(x_grid, x_val_clamped)
    j = searchsortedlast(t_grid, t_val_clamped)
    
    # Handle boundaries
    i = min(i, length(x_grid) - 1)
    j = min(j, length(t_grid) - 1)
    
    # Get the four corner values
    v11 = values[i, j]
    v12 = values[i, j+1]
    v21 = values[i+1, j]
    v22 = values[i+1, j+1]
    
    # Get grid points
    x1, x2 = x_grid[i], x_grid[i+1]
    t1, t2 = t_grid[j], t_grid[j+1]
    
    # Bilinear interpolation:
    # 1. Interpolate in x direction for both t boundaries
    # 2. Interpolate in t direction between the two x-interpolated values
    
    # For extrapolation, use clamped values for weight calculation
    # but preserve Dual number structure by using original x, t when within bounds
    # When clamped, we need to use the clamped value for proper extrapolation
    if x_val == x_val_clamped
        x_for_interp = x
    else
        # Convert clamped value to same type as x to preserve Dual structure
        x_for_interp = convert(typeof(x), x_val_clamped)
    end
    
    if t_val == t_val_clamped
        t_for_interp = t
    else
        # Convert clamped value to same type as t to preserve Dual structure
        t_for_interp = convert(typeof(t), t_val_clamped)
    end
    
    # Compute weights (preserves Dual number structure)
    wx = (x_for_interp - x1) / (x2 - x1)
    wt = (t_for_interp - t1) / (t2 - t1)
    
    # Interpolate in x direction at t1 and t2
    v_t1 = v11 + wx * (v21 - v11)  # value at (x, t1)
    v_t2 = v12 + wx * (v22 - v12)  # value at (x, t2)
    
    # Interpolate in t direction
    return v_t1 + wt * (v_t2 - v_t1)
end

"""
    Flat()

Placeholder type for extrapolation boundary condition (flat extrapolation).
Used for API compatibility with Interpolations.jl.
"""
struct Flat end

"""
    linear_interpolation(grid, values; extrapolation_bc=Flat())

Create a simple linear interpolation object that is fully compatible with ForwardDiff.jl.

# Arguments
* `grid`: Tuple of grid points. For 1D: `(grid_points,)`, for 2D: `(x_grid, t_grid)`
* `values`: Array of values at grid points. Vector for 1D, Matrix for 2D.
* `extrapolation_bc`: Extrapolation boundary condition (only `Flat()` is supported)

# Returns
A `SimpleLinearInterpolation` object that can be called like a function.

# Examples
```julia
# 1D interpolation
itp1d = linear_interpolation(([0.0, 1.0, 2.0],), [10.0, 20.0, 30.0])
value = itp1d(0.5)  # Returns 15.0

# 2D interpolation
x_grid = [0.0, 1.0, 2.0]
t_grid = [0.0, 10.0, 20.0]
values = [10.0 20.0 30.0; 15.0 25.0 35.0; 20.0 30.0 40.0]
itp2d = linear_interpolation((x_grid, t_grid), values)
value = itp2d(0.5, 5.0)  # Returns interpolated value
```
"""
function linear_interpolation(grid::Tuple, values::AbstractArray; extrapolation_bc=Flat())
    D = length(grid)
    
    if D == 1
        # 1D interpolation
        grid_vec = collect(grid[1])
        values_vec = collect(values)
        T = promote_type(eltype(grid_vec), eltype(values_vec))
        return SimpleLinearInterpolation{1, T}((convert(Vector{T}, grid_vec),), convert(Vector{T}, values_vec))
    elseif D == 2
        # 2D interpolation
        grid1 = collect(grid[1])
        grid2 = collect(grid[2])
        T = promote_type(eltype(grid1), eltype(grid2), eltype(values))
        return SimpleLinearInterpolation{2, T}(
            (convert(Vector{T}, grid1), convert(Vector{T}, grid2)),
            convert(Array{T, 2}, values)
        )
    else
        error("Only 1D and 2D interpolation are supported")
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

#---End-Simple-Linear-Interpolation-Implementation---