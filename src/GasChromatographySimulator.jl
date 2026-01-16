module GasChromatographySimulator

using Reexport
using Integrals
@reexport using OrdinaryDiffEq: OwrenZen3, OwrenZen4, OwrenZen5, ODEProblem
using ForwardDiff
using ForwardDiffOverMeasurements
using DataFrames
using CSV
using Plots
using HypertextLiteral
using PlutoUI
using ChemicalIdentifiers
using UrlDownload
using Measurements

# some constants
const Tst = 273.15            # K
const R = 8.31446261815324    # J mol⁻¹ K⁻¹
const Tn = 25.0 + Tst         # K
const pn = 101300             # Pa
const custom_database_url = "https://raw.githubusercontent.com/JanLeppert/RetentionData/main/data/add_CI_db.tsv"
const shortnames_url = "https://raw.githubusercontent.com/JanLeppert/RetentionData/main/data/shortnames.csv"
const missing_url = "https://raw.githubusercontent.com/JanLeppert/RetentionData/main/data/missing.csv"
const shortnames_filepath = joinpath("data", "shortnames.csv")
const missing_filepath = joinpath("data", "missing.csv")
#const k_th = 1e20#1e12            # threshold for retention factor k, if k>k_th => k=k_th 

# Loading files
## Structures and Constructors
include("./Structures.jl")
## Physical-model-functions
include("./Model.jl")
## solving-functions
include("./Solving.jl")
## UI-functions
include("./UI.jl")
## plot-functions
include("./Plot.jl")
## misc-functions
include("./Misc.jl")



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
    
    # Compute weights (preserves Dual number structure)
    wx = (x - x1) / (x2 - x1)
    wt = (t - t1) / (t2 - t1)
    
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

#---End-Simple-Linear-Interpolation-Implementation---

#---Begin-Functions-used-for-Parameter-construction 
"""
    gradient(x, a; Tcontrol="inlet")

Defines the gradient of the column diameter, film thickness or
temperature along the GC column.  

# Arguments
* `x`: Position along the GC column, in m.
* `a`: Parameters of the gradient function

The form of `a` decides the actual used function for the gradient. The
following options are realized:
* `a` is a single value (e.g. Float or Int): The gradient function is constant for all positions `x` with the value of `a`. 
* `a` is a 1D-array of length = 4: The gradient function is a exponential function. The 4 values of `a`: `a = [f₀, x₀, L₀, α]` with 
    * `f₀`: Start value at `x = x₀`.
    * `x₀`: Shift in `x` position.
    * `L₀`: Distance over which the value drops from f₀ to 0.
    * `α`: Factor describing the gradient profile.
    if `α <= 0` 
        `f = f₀ .* (1 .- exp.(α.*(1 .- (x.-x₀) ./ L₀)) + (1 .- (x.-x₀) ./ L₀)
        .* exp.(α))`
    if `α > 0`
        `f = f₀ .* (exp.(-α.*(x.-x₀) ./ L₀) .- (x.-x₀) ./ L₀ .* exp.(-α))` 
* `a` is a 2D-array with size = (n, 4): The gradient function is an
    exponential function. The 4 entrys have the same meaning as above, but their
    values can change over the times defined in the `time_steps` of the
    Program structure. At these `time_steps[i]` the gradient function is
    described by the corresponding parameters `a[i,:]`. Between the
    `time_steps[i]` and `time_steps[j]` the value of the gradient function at position `x` is linear
    interpolated from the gradien functions defined by `a[i,:]` and `a[j,:]`.

# Examples
```julia
julia> d(x) = gradient(x, [0.1e-3])
```

```julia
julia> gf(x) = gradient(x, [[20.0, 20.0, 20.0, 20.0] zeros(4) 10.0.*ones(4) [0.0, -2.0, -5.0, -5.0]])
```
"""
# gradient functions
function gradient(x, a; Tcontrol="inlet")
    if size(a) == (1,)
        # constant value, no gradient
        f = a[1]
    elseif size(a) > (1,)
        if  length(size(a))==1
            if length(a)==4
                # for diameter or film thickness, values of parameters 'a' are fixed
                # over time
                f₀ = a[1] # start value
                x₀ = a[2] # shift in x 
                L₀ = a[3] # length of the linear segment
                α  = a[4] # exponential factor, α=0 -> linear
                if α<=0 # concave/linear, weak change at beginning and strong change at end of column
                    f = f₀ .* (1 .- exp.(α.*(1 .- (x.-x₀) ./ L₀)) + (1 .- (x.-x₀) ./ L₀) .* exp.(α))
                elseif α>0 # convex function, strong change at beginning and weak change at end of the column
                    f = f₀ .* (exp.(-α.*(x.-x₀) ./ L₀) .- (x.-x₀) ./ L₀ .* exp.(-α))
                end
            # other functions ...
            end
        elseif length(size(a)) == 2
            if a==zeros(size(a))
                f = zeros(size(a)[1])
                #if size(a)[2] == 1
                # for thermal gradient, no change in time  of values of the
                # parameters
            elseif size(a)[2] == 4
                # for thermal gradient, values of parameters 'a[i,:]' can change
                # over time
                f₀ = a[:,1] # start value
                x₀ = a[:,2] # shift in x (e.g. to correct for the real position of IR-sensors)
                L₀ = a[:,3] # length of the linear segment
                α  = a[:,4] # exponential factor, α=0 -> linear
                f = Array{Real}(undef, length(α))
                for i=1:length(α)
                    if α[i]<=0 # concave function, weak change at beginning and strong change at end of column
                        if Tcontrol=="inlet"
                            f[i] = f₀[i] .* (.- exp.(α[i].*(1 .- (x.-x₀[i]) ./ L₀[i])) + (1 .- (x.-x₀[i]) ./ L₀[i]) .* exp.(α[i]))
                        elseif Tcontrol=="outlet"
                            f[i] = f₀[i] .* (1 .- exp.(α[i].*(1 .- (x.-x₀[i]) ./ L₀[i])) + (1 .- (x.-x₀[i]) ./ L₀[i]) .* exp.(α[i]))
                        end
                    elseif α[i]>0 # convex function, strong change at beginning and weak change at end of the column
                        if Tcontrol=="inlet"
                            f[i] = f₀[i] .* (exp.(-α[i].*(x.-x₀[i]) ./ L₀[i]) .- (x.-x₀[i]) ./ L₀[i] .* exp.(-α[i]) .- 1)
                        elseif Tcontrol=="outlet"
                            f[i] = f₀[i] .* (exp.(-α[i].*(x.-x₀[i]) ./ L₀[i]) .- (x.-x₀[i]) ./ L₀[i] .* exp.(-α[i]))
                        end
                    end
                end
            # other functions ...
            end
        end
    end
    return f
end

"""
    temperature_interpolation(time_steps, temp_steps, gradient_function, L)

Construct the temperature function depending on position `x` and
time `t`.  

# Arguments
* `time_steps::Array{<:Real,1}`: Time steps in s (seconds). 
* `temp_steps::Array{<:Real,1}`: Temperature steps in °C (degree celsius).
* `gf::Function`: Gradient function `gf(x, a_gf)`, describes the thermal gradient.
* `L`: Length of the capillary measured in m (meter).

For the spatial dependency of the interpolated temperature `T_ipt(x,t)` the
gradient function `gf` is calculated every 1e-3 m (1 mm). Positions
inbetween are linear interpolated. For the temporal dependency the
temperatures `temp_steps` defined at the `time_steps` are linear
interpolated over time `t`.   

# Examples
```julia
julia> T_itp = temperature_interpolation([0.0, 60.0, 300.0, 120.0], [40.0, 40.0, 320.0, 320.0], gf, 10.0)
```
"""
function temperature_interpolation(time_steps::Array{<:Real,1}, temp_steps::Array{<:Real,1}, gradient_function::Function, L; ng=false, dx=1e-3)
	T(x) = temp_steps .+ gradient_function(x) 
	L_ = Measurements.value(L)
    if ng == false
	    nx = 0.0:dx:L_ # mm exact
    else
        nx = [0.0, L_]
    end
	nt = cumsum(time_steps)
	Tmat = Array{Real}(undef, length(nx), length(nt))
	for j=1:length(nt)
		for i=1:length(nx)
			Tmat[i,j] = T(nx[i])[j] + 273.15
		end
	end
	T_itp = linear_interpolation((collect(nx), nt), Tmat, extrapolation_bc=Flat())
	return T_itp
end

"""
    steps_interpolation(time_steps, steps)

Construct a linear interpolated function depending on time `t` of the `steps`-values over `time_steps`.  

# Arguments
* `time_steps::Array{<:Real,1}`: Time steps in s (seconds). 
* `steps::Array{<:Real,1}`: steps, e.g. pressure or flow. 

# Examples
```julia
julia> pin_itp = steps_interpolation([0.0, 60.0, 300.0, 120.0], 
                                    [300000.0, 300000.0, 400000.0, 400000.0])
```
"""
function steps_interpolation(time_steps::Array{<:Real,1}, steps::Array{<:Real,1})
    s_itp = linear_interpolation((cumsum(time_steps), ), steps, extrapolation_bc=Flat())
    return s_itp
end

"""
    CAS_identification_from_CAS(cas_number::AbstractString)

Look up a substance by its CAS number using ChemicalIdentifiers.jl to find the name, formula, molecular weight `MW` and `smiles`-identifier. 
Returns a NamedTuple with the substance information or placeholder data if the CAS number is not found.

# Example
```julia
julia> CAS_identification_from_CAS("71-43-2")
(Name = "Benzene", CAS = "71-43-2", formula = "C6H6", MW = 78.11, smiles = "c1ccccc1")
```
"""
function CAS_identification_from_CAS(cas_number::AbstractString; show_info=false)
    # Check if the input matches CAS number format (XX-XX-X)
    if !occursin(r"^\d{1,7}-\d{2}-\d$", cas_number)
        error("Invalid CAS number format. Expected format: XX-XX-X")
    end

    # Split the CAS number into its components
    cas_parts = split(cas_number, "-")
    cas_tuple = (parse(Int, cas_parts[1]), parse(Int, cas_parts[2]), parse(Int, cas_parts[3]))

    if show_info
        @info "Searching for CAS number $cas_number"
    end
    ci = try
        search_chemical(cas_tuple)
    catch
        missing
    end
    if show_info
        @info "Found $ci"
    end

    if ismissing(ci)
        # Placeholder data if CAS number not found
        id = (Name = string(cas_number,"_ph"), 
                CAS = "629-62-9_ph", 
                formula = "C15H32", 
                MW = 212.41, 
                smiles = "CCCCCCCCCCCCCCC")
    else
        # Reconstruct CAS number with proper formatting
        if length(digits(ci.CAS[2])) == 1
            CAS = string(ci.CAS[1], "-0", ci.CAS[2], "-", ci.CAS[3])
        else
            CAS = string(ci.CAS[1], "-", ci.CAS[2], "-", ci.CAS[3])
        end
        id = (Name = ci.iupac_name, 
              CAS = CAS, 
              formula = ci.formula, 
              MW = ci.MW, 
              smiles = ci.smiles)
    end
    if show_info
        @info "Found $id"
    end
    return id
end

"""
	CAS_identification(Name)

Look up the substance name from the `data` dataframe with ChemicalIdentifiers.jl to find the `CAS`-number, the `formula`, the molecular weight `MW` and the `smiles`-identifier. If the name is not found in the database of ChemicalIdentifiers.jl a list with alternative names (`shortnames.csv`) is used. If there are still no matches, `missing` is used.
"""
function CAS_identification(Name; show_info=false)
    # Check if the input is in CAS number format
    if occursin(r"^\d{1,7}-\d{2}-\d$", Name)
        return CAS_identification_from_CAS(Name)
    else
        load_custom_CI_database(custom_database_url)
        shortnames = try
            DataFrame(urldownload(shortnames_url))
        catch
            DataFrame(CSV.File(shortnames_filepath))
        end
        missingsubs = try 
            DataFrame(urldownload(missing_url))
        catch
            DataFrame(CSV.File(missing_filepath))
        end
        if Name in missingsubs.name # first look the name up in the missing.csv
            j = findfirst(Name.==missingsubs.name)
            id = (Name = missingsubs.name[j], CAS = missingsubs.CAS[j], formula = missingsubs.formula[j], MW = missingsubs.MW[j], smiles = missingsubs.smiles[j])
        else
            if Name in shortnames.shortname # second look for name alternative in shortnames.csv for an alternative name used for the search
                j = findfirst(Name.==shortnames.shortname)
                name = String(shortnames.name[j])
            else # otherwise use the input name
                name = Name
            end
            if show_info
                @info "Searching for $name"
            end
            ci = try
                search_chemical(name)
            catch
                missing
            end
            if show_info
                @info "Found $ci"
            end
            if ismissing(ci) # 
                #ci_ = search_chemical("629-62-9")
                # Placeholder data from Pentadecane
                id = (Name = string(Name,"_ph"), 
                        CAS = "629-62-9_ph", 
                        formula = "C15H32", 
                        MW = 212.41, 
                        smiles = "CCCCCCCCCCCCCCC")
            else
                if length(digits(ci.CAS[2])) == 1 # if the second CAS number has only one digit, add a leading zero
                    CAS = string(ci.CAS[1], "-0", ci.CAS[2], "-", ci.CAS[3])
                else
                    CAS = string(ci.CAS[1], "-", ci.CAS[2], "-", ci.CAS[3])
                end
                id = (Name = Name, CAS = CAS, formula = ci.formula, MW = ci.MW, smiles = ci.smiles)
            end
        end
        if show_info
            @info "Found $id"
        end
        return id
    end
end

"""
    load_solute_database(db_path, db, sp, gas, solutes, t₀, τ₀) 

Load the data of `solutes` for the stationary phase `sp` and the mobile
phase `gas` from the database file `db` (located in `db_path`) into an array
of the structure `Substance`. The row number of the selected solutes in the loaded database are added to the annotations of the `Substance` structure. 

# Arguments
* `db_path::String`: Path to the database file.
* `db::String`: Name of the database file. 
* `sp::String`: Name of the stationary phase.
* `gas::String`: Name of the mobile phase.
* `solutes::Array{<:AbstractString,1}`: Name of the solutes.
* `t₀::Array{<:Real,1}`: Initial start times of the solutes.
* `τ₀::Array{<:Real,1}`: Initial peak widths of the solutes. 

# Examples
```julia
julia> sub = load_solute_database("path/to/the/file", "db.csv", "DB5", "He", ["C10", "C11"], [0.0, 0.0], [0.5, 0.5])
```
"""
function load_solute_database(db_path::String, db::String, sp::String, gas::String, solutes::Array{<:AbstractString,1}, t₀::Array{<:Real,1}, τ₀::Array{<:Real,1})
	# load the information about a solute from a data base and collecting these informations in the 
    # structure Substance
    # 
	# db_path ... Path to the database file
	# db ... Name of the database
	# sp ... Name of the stationary phase
	# gas ... Name of the used mobile phase gas
	# solutes ... Names of the solutes for which the informations should be used
    # τ₀ ... initial values of the peak width for the solutes
    # t₀ ... initial values of the time for the solutes
	db_dataframe = DataFrame(CSV.File(string(db_path,"/",db), header=1, silencewarnings=true, stringtype=String))
    insertcols!(db_dataframe, 1, :No => collect(1:length(db_dataframe.Name)))
    sub = load_solute_database(db_dataframe, sp, gas, solutes, t₀, τ₀)
	return sub
end

"""
	find_index_in_database(solutes, db)

Find the index of the substances `solutes` in the database `db`. Herby `solutes` can be a 
* vector of integers with the number `No` of the substances in the database.
* vector of strings containing substance names.
* vector of strings in the CAS number format.

The index of the found substance in the database `db` and the index of the not found substances in the `solutes` vector will be returned.
"""
function find_index_in_database(solutes, db)
	# find the index in the database
	# three cases
	# 1. number No
	# 2. CAS number
	# 3. name
	regex_CAS = r"([0-9]{2,8}-[0-9]{2}-[0-9]{1})"
	if isa(solutes, Array{Int, 1}) # number No
		colname = :No
	elseif isa(solutes, Array{String, 1}) #  name
		if all(occursin.(regex_CAS, solutes)) # CAS number
			colname = :CAS
		else
			colname = :CAS # if name in solutes is different from db it will not find it -> go with CAS number
			id = GasChromatographySimulator.CAS_identification.(solutes)
			solutes = [id[x].CAS for x in 1:length(id)]
		end
	else
		error("selected substances should be either a vector of integers (number No. of the substance in the database), a vector of strings in the CAS-number format or a vector of strings containing substance names.")
	end
	# all indices of the solutes in the database:
	index = [findall(solutes[x].==db[!, colname]) for x in 1:length(solutes)]
	# transform to vector:
	i_db = filter(isinteger, reduce(vcat, index))
	# indices of the solutes in the solutes-vector, which are not found in the database:
	i_solutes_not_found = findall(isempty.(index))
	# indices of the solutes in the solutes-vector of the found database indices:
	i_db_i_solutes = reduce(vcat,[findall(db[!, colname][reduce(vcat, index)][x].==solutes) for x in 1:length(reduce(vcat, index))])
return i_db, i_solutes_not_found, i_db_i_solutes
end

"""
	annotations(db)

Extracts an annotation for every substance from the database `db`. The annotation consists of the entries for:
* `Source`
* optional categories `Cat` (all columns with the name starting with "Cat")
* optional number `No`
"""
function annotations(db)
	i_Catnames = findall(occursin.("Cat", names(db)))
    nCat = length(i_Catnames)
    if nCat < 1
        Annotation = String.(db.Source)
    else
        Annotation = String.(db.Source)
        for i=1:nCat
        	for j=1:length(Annotation)
                if typeof(db[j,i_Catnames[i]]) != Missing
                    Annotation[j] = string(Annotation[j], ", ", db[j,i_Catnames[i]])
                end
            end
        end
    end
	if "No" in names(db)
		Annotation = Annotation.*string.(", No: ", db.No)
	end
	return Annotation
end

# new version
"""
	load_solute_database(db_, sp_, gas_, solutes_, t₀_, τ₀_)

New function to load the data for the substances `solutes_` from the database `db_` for the stationary phase `sp_` and mobile phase `gas_` with initial time `t₀_` and initial peak width `τ₀_`. The parameters `solutes_`, `t₀_` and `τ₀_` must be vectors of the same length. Acceptable values for `solutes_` are either the substance names, CAS numbers or the number No from the database.

The loaded data is returned as a vector of the GasChromatographySimulator.Substance structure.
"""
function load_solute_database(db_, sp_, gas_, solutes_, t₀_, τ₀_)
	# check for same size of solutes_, t₀_ and τ₀_:
	if (length(solutes_) != length(t₀_)) || (length(solutes_) != length(τ₀_))
		@warn "Number of solutes ($(length(solutes_))), t₀ ($(length(t₀_))) and τ₀ ($(length(τ₀_))) do not match. Missing values of t₀ or τ₀ will be added with value 0.0 s. Additional values of t₀ or τ₀ will be skipped."
        t₀_ = same_length(t₀_, solutes_)
        τ₀_ = same_length(τ₀_, solutes_)
	end	
    ## remove solutes with missing CAS and give warning
    if true in ismissing.(db_.CAS) 
        @warn "Some CAS-numbers are missing. These substances are skipped."
        db = disallowmissing!(db_[completecases(db_, :CAS), :], :CAS)
    else
        db = db_
    end # is this still needed???

	if sp_ == "" # no stationary phase is selected, like for transferlines
		i_db, i_solutes_not_found, i_db_i_solutes = find_index_in_database(solutes_, db)
		if isempty(i_solutes_not_found) == false
			@warn "Some solutes could not be found: $(solutes_[i_solutes_not_found])."
		end
       	Name = db.Name[i_db]
        CAS = db.CAS[i_db]
        # use placeholder values
        Tchar = ones(length(Name))
        θchar = ones(length(Name))
        ΔCp = ones(length(Name))
        φ₀ = ones(length(Name))
        Annotation = fill("no sp", length(Name))
		# Cag -> reading or calculating
		i_Cag = findfirst("Cag".==names(db))
        # parameter name + "_" indicates the column with uncertainty values of this quantity 
		i_Cag_unc = findfirst(occursin.("Cag_", names(db)))
		if isa(i_Cag, Int)
			if isa(i_Cag_unc, Int)
				Cag = db[i_db, i_Cag] .± db[i_db, i_Cag_unc]
			else
				Cag = db[i_db, i_Cag]
			end
		else
			id = GasChromatographySimulator.CAS_identification.(CAS)
			Cag = GasChromatographySimulator.diffusivity.(id, gas_)
		end
       	t₀ = t₀_[i_db_i_solutes]
        τ₀ = τ₀_[i_db_i_solutes]
	else
		# 1. Filter the stationary phase
        db_filtered = filter([:Phase] => x -> x==sp_, db)
		if isempty(db_filtered)
			@warn "No data for selected solutes $(solutes_) for stationary phase $(sp_)."
			#extracted_db = DataFrame()
			i_db = Int[]
			i_solutes_not_found = collect(1:length(solutes_))
			sub = GasChromatographySimulator.Substance[]
		else
			i_db, i_solutes_not_found, i_db_i_solutes = find_index_in_database(solutes_, db_filtered)
			if isempty(i_solutes_not_found) == false
				@warn "Some solutes could not be found for stationary phase $(sp_): $(solutes_[i_solutes_not_found])."
			end
			
			# look for the following (optional) columns in the data 
			# parameter name + "_" indicates the column with uncertainty values of this quantity 
			i_Tchar_unc = findfirst(occursin.("Tchar_", names(db_filtered)))
			i_θchar_unc = findfirst(occursin.("thetachar_", names(db_filtered)))
			i_ΔCp = findfirst("DeltaCp".==names(db_filtered))
			i_ΔCp_unc = findfirst(occursin.("DeltaCp_", names(db_filtered)))
			i_Cag = findfirst("Cag".==names(db_filtered))
			i_Cag_unc = findfirst(occursin.("Cag_", names(db_filtered)))
			i_φ₀_unc = findfirst(occursin.("phi0_", names(db_filtered)))
	
			Name = db_filtered.Name[i_db]
	        CAS = db_filtered.CAS[i_db]
			# Tchar
			Tchar = if isa(i_Tchar_unc, Int)
				db_filtered.Tchar[i_db] .± db_filtered[i_db, i_Tchar_unc] .+ GasChromatographySimulator.Tst
			else
				db_filtered.Tchar[i_db] .+ GasChromatographySimulator.Tst
			end
			# θchar
			θchar = if isa(i_θchar_unc, Int)
				db_filtered.thetachar[i_db] .± db_filtered[i_db, i_θchar_unc]
			else
				db_filtered.thetachar[i_db]
			end
			# ΔCp
	        ΔCp = if isa(i_ΔCp, Int)
				if isa(i_ΔCp_unc, Int)
					db_filtered[i_db, i_ΔCp] .± db_filtered[i_db, i_ΔCp_unc]
				else
					db_filtered[i_db, i_ΔCp]
				end
			else
				ones(length(i_db))
			end
			# phi0
			φ₀ = if isa(i_φ₀_unc, Int)
				db_filtered.phi0[i_db] .± db_filtered[i_db, i_φ₀_unc]
			else
				db_filtered.phi0[i_db]
			end
			# annotations
	        Annotation = annotations(db[i_db,:]) 
			# Cag -> reading or calculating
			Cag = if isa(i_Cag, Int)
				if isa(i_Cag_unc, Int)
					db_filtered[i_db, i_Cag] .± db_filtered[i_db, i_Cag_unc]
				else
					db_filtered[i_db, i_Cag]
				end
			else
				GasChromatographySimulator.diffusivity.(GasChromatographySimulator.CAS_identification.(CAS), gas_)
			end
	        t₀ = [t₀_[i_db_i_solutes[x]] for x in 1:length(i_db_i_solutes)]
	        τ₀ = [τ₀_[i_db_i_solutes[x]] for x in 1:length(i_db_i_solutes)]
		end
	end
	sub = [GasChromatographySimulator.Substance(Name[x],
                            CAS[x],
                            Tchar[x], 
                            θchar[x], 
                            ΔCp[x], 
                            φ₀[x],
                            Annotation[x],
                            Cag[x], 
                            t₀[x],
                            τ₀[x]) for x in 1:length(i_db)]
	return unique(sub) # remove duplicates
end

"""
    all_solutes(sp, db; id=false, separator=" - ") 

Extract the name of all solutes for which data in a database `db` and the
stationay phase `sp` is available. 

# Arguments
* `sp`: Name of the stationary phase.
* `db`: DataFrame of the database.
* `id`: if `true` than the number of the solutes in the database are combined with the solute name. Default = false.
* `separator`: string used to separate the solute number from the solute name. Default = " - ". 

# Examples
```julia
julia> all = all_solutes("DB5", db)
```
"""
function all_solutes(sp, db; id=false, sperator=" - ")
	db_filter = filter([:Phase] => x -> x==sp, db)
	if id == true
        if "No" in names(db)
		    solutes = string.(db_filter.No, sperator, db_filter.Name)
        else
            @warn "For 'id=true' a column 'No' with the id numbers is needed in the database 'db'."
        end
	else
		solutes = string.(db_filter.Name)
	end
	return solutes
end

"""
    conventional_program(CP; time_unit="min")

Translate the conventional program notation into a vector of time steps and value steps (temperature, pressure, flow) used in GasChromatographySimulator.Program

The conventional temperature program is defined as an array of the following form (for a temperature program):
`CP = [T₁, t₁, r₁, T₂, t₂, r₂, T₃, t₃, ...]` corresponding to the notation:
`T₁(t₁) - r₁ - T₂(t₂) - r₂ - T₃(t₃) - ...` which can be read as:
Starting temperature `T₁` is holded for time `t₁`. After the holding time the temperature increases linearly with the heating rate `r₁`, until temperature `T₂` is reached. This temperature is held for the time `t₂` after which the temperature increases linearly by the heating rate `r₂` until temperature `T₃` is reached, which is hold for the time `t₃`, and so on. 

The option `time_unit` determines the unit of time in the program `CP`. For `time_unit = "min"` (default) the times are measured in minutes and the heating rates in °C/min. For `time_unit = "s"` the times are measured in seconds and the heating rate in °C/s. 
"""
function conventional_program(CP; time_unit="min")
    if time_unit == "min"
        c = 60.0
    elseif time_unit == "s"
        c = 1.0
    end
    value_steps = Real[]
    for i=1:3:length(CP) # every third CP-entry starting from 1 is a value step
        push!(value_steps, CP[i])
        push!(value_steps, CP[i])
    end
    hold_times = Real[]
    for i=2:3:length(CP)
        push!(hold_times, CP[i])
    end
    heating_rates = Real[]
    for i=3:3:length(CP)
        push!(heating_rates, CP[i])
    end
    time_steps = Array{Real}(undef, length(value_steps))
    time_steps[1] = 0.0
    for i=2:2:length(value_steps)
        time_steps[i] = hold_times[Int(i/2)] * c
        if i<length(value_steps)
            time_steps[i+1] = (value_steps[i+1] - value_steps[i])/heating_rates[Int(i/2)] * c
        end
    end
    if time_steps[1] == time_steps[2] && value_steps[1] == value_steps[2]
		time_steps = time_steps[2:end]
		value_steps = value_steps[2:end]
	end
    return time_steps, value_steps
end

"""
    temperature_program(time_steps, value_steps; time_unit="min")

Translate the vector of time steps and value steps (temperature, pressure, flow) into a conventional program notation.

The conventional temperature program is defined as an array of the following form (for a temperature program):
`CP = [T₁, t₁, r₁, T₂, t₂, r₂, T₃, t₃, ...]` corresponding to the notation:
`T₁(t₁) - r₁ - T₂(t₂) - r₂ - T₃(t₃) - ...` which can be read as:
Starting temperature `T₁` is holded for time `t₁`. After the holding time the temperature increases linearly with the heating rate `r₁`, until temperature `T₂` is reached. This temperature is held for the time `t₂` after which the temperature increases linearly by the heating rate `r₂` until temperature `T₃` is reached, which is hold for the time `t₃`, and so on. 

The option `time_unit` determines the unit of time in the program `CP`. For `time_unit = "min"` (default) the times are measured in minutes and the heating rates in °C/min. For `time_unit = "s"` the times are measured in seconds and the heating rate in °C/s. 
"""
function temperature_program(time_steps, value_steps; time_unit="min")
    if time_unit == "min"
        c = 60.0
    elseif time_unit == "s"
        c = 1.0
    end
    # identify temperature pairs (the same value of temperature at neigboring elements)
    index_pair = Int[]
    for i=1:(length(value_steps)-1)
        if value_steps[i] == value_steps[i+1]
            push!(index_pair, i)
        end
    end
    # identify single temperatures
    index_single = Int[]
    for i=1:length(value_steps)
        if (i in index_pair) == false && (i in index_pair.+1) == false
            push!(index_single, i)
        end
    end
    values = value_steps[sort([index_pair; index_single])]
    # every (1+(i-1)*3)th element of VP is a value element of `values`
    
    thold = Array{Real}(undef, length(values))
    a = sort([index_pair.+1; index_single])
    for i=1:length(values)
        if a[i] in index_single # holding time is zero for single entrys
            thold[i] = 0.0
        else # holding times are for paired temperatures the time_steps[index_pair.+1]
            thold[i] = time_steps[a[i]] / c
        end
    end 
    # every (2+(i-1)*3)th element of VP is a holding time

    rate = Array{Real}(undef, length(values)-1)
    theat = time_steps[sort([index_pair; index_single])[2:end]]
    for i=1:(length(values)-1)
        rate[i] = (values[i+1] - values[i])/theat[i] * c
    end

    VP = Array{Real}(undef, 2 + (length(values)-1)*3)
    for i=1:length(values)
        VP[1+(i-1)*3] = values[i]
    end
    for i=1:length(thold)
        VP[2+(i-1)*3] = thold[i]
    end
    for i=1:length(rate)
        VP[3+(i-1)*3] = rate[i]
    end
    return VP
end

"""
    common_time_steps(time_steps_1, time_steps_2)

Estimate a new set of time steps, which represents the combination of `time_steps_1` and `time_steps_2`.
"""
function common_time_steps(time_steps_1, time_steps_2)
	# constructs the new timesteps common for all moduls
	ctselements = Real[]
    for j=1:length(time_steps_1)
        push!(ctselements, cumsum(time_steps_1)[j])
    end
    for j=1:length(time_steps_2)
        push!(ctselements, cumsum(time_steps_2)[j])
    end
	new_time_steps = [0.0; diff(sort(unique(ctselements)))]
	return new_time_steps
end

"""
    new_value_steps(value_steps, time_steps, new_time_steps)

Estimate the new value steps at the `new_time_steps` from the original set of `value_steps` over `time_steps`. The new values at new time steps are calculated from a linear change of the value between the original time steps.
""" 
function new_value_steps(value_steps, time_steps, new_time_steps)
    new_values = Array{Real}(undef, length(new_time_steps))
    index_calc = Int[]
    index_old = Int[]
    for i=1:length(new_time_steps)
        if cumsum(new_time_steps)[i] in cumsum(time_steps)
            j = findfirst(cumsum(new_time_steps)[i].==cumsum(time_steps))
            new_values[i] = value_steps[j]
            push!(index_old, i)
        else
            push!(index_calc, i)
        end
    end
    for i=1:length(index_calc)
        i1 = findlast(index_old.<index_calc[i])
        i2 = findfirst(index_old.>index_calc[i])
        if isnothing(i1) || isnothing(i2)
            new_values[index_calc[i]] = value_steps[end]
        else
            if time_steps[i2] == 0.0
                i2 = i2 + 1
            end
            v1 = value_steps[i1]
            v2 = value_steps[i2]
            t1 = cumsum(time_steps)[i1]
            t2 = cumsum(time_steps)[i2]
            rate = (v2 - v1)/(t2 - t1)
            if t1 == t2
                new_values[index_calc[i]] = v2
            else
                new_values[index_calc[i]] = v1 + rate * (cumsum(new_time_steps)[index_calc[i]] - t1)
            end
        end
    end
    return new_values
end
#---End-Functions-used-for-Parameter-construction--- 

#---Begin-Result-Functions---
"""
    peaklist(sol, par)

Construct a DataFrame with the peak list of the solution `sol` of the
simulation of the GC system defined by `par`. 

# Output

The peaklist DataFrame consists of the entrys:
* `No`: Number of the solute in the database. 
* `Name`: Name of the solute.
* `tR`: Retention time of the solute (in s).
* `τR`: Peak width of the solute (in s). 
* `TR`: Temperature of the end of the column at the retention time (in °C).
* `σR`: Band width of the solute at retention time (in m).
* `uR`: Solute velocity at retention time (in m/s).
* `kR`: Retention factor of the solute at retention time.
* `Res`: Resolution (4τ) between neighboring peaks.
* `Δs`: separation metric between neighboring peaks, assuming linear development of peak width `τR` between the peaks.
* `Annotations`: additional anotations, e.g. Source, categories, if available

# Examples

```julia
julia> pl = peaklist(sol, par)
...
```    
"""
function peaklist(sol, par; thread=false)
    if thread == true
        pl = peaklist_thread(sol, par)
    else
        pl = peaklist_unthread(sol, par)
    end
    return pl
end

function peaklist(sol, L, d, df, gas, T_itp, Fpin_itp, pout_itp, Name, CAS, ann, Tchar, θchar, ΔCp, φ₀, opt; thread=false)
    if thread == true
        pl = peaklist_thread(sol, L, d, df, gas, T_itp, Fpin_itp, pout_itp, Name, CAS, ann, Tchar, θchar, ΔCp, φ₀, opt)
    else
        pl = peaklist_unthread(sol, L, d, df, gas, T_itp, Fpin_itp, pout_itp, Name, CAS, ann, Tchar, θchar, ΔCp, φ₀, opt)
    end
    return pl
end

function peaklist_thread(sol, par)
	n = length(par.sub)
    T = typeof(sol[1].u[end][1])
    # sol is solution from ODE system
    No = Array{Union{Missing, Int64}}(undef, n)
    Name = Array{String}(undef, n)
    CAS = Array{String}(undef, n)
    tR = Array{T}(undef, n)
    TR = Array{T}(undef, n)
    σR = Array{T}(undef, n)
    uR = Array{T}(undef, n)
    τR = Array{T}(undef, n)
    kR = Array{T}(undef, n)
    Res = fill(NaN, n)
    Δs = fill(NaN, n)
    Annotations = Array{String}(undef, n)
    Threads.@threads for i=1:n
        Name[i] = par.sub[i].name
        CAS[i] = par.sub[i].CAS
        if sol[i].t[end]==par.col.L
            tR[i] = sol[i].u[end][1]
            TR[i] = par.prog.T_itp(par.col.L, tR[i]) - 273.15 
            uR[i] = 1/residency(par.col.L, tR[i], par.prog.T_itp, par.prog.Fpin_itp, par.prog.pout_itp, par.col.L, par.col.d, par.col.df, par.col.gas, par.sub[i].Tchar, par.sub[i].θchar, par.sub[i].ΔCp, par.sub[i].φ₀; ng=par.opt.ng, vis=par.opt.vis, control=par.opt.control, k_th=par.opt.k_th)
            τR[i] = sqrt(sol[i].u[end][2])
            σR[i] = τR[i]*uR[i]
            kR[i] = retention_factor(par.col.L, tR[i], par.prog.T_itp, par.col.d, par.col.df, par.sub[i].Tchar, par.sub[i].θchar, par.sub[i].ΔCp, par.sub[i].φ₀; k_th=par.opt.k_th)
        else
            tR[i] = NaN
            TR[i] = NaN
            uR[i] = NaN
            τR[i] = NaN
            σR[i] = NaN
            kR[i] = NaN
        end
        No[i] = try
            parse(Int,split(par.sub[i].ann, ", ")[end])
        catch
            try 
                parse(Int,split(par.sub[i].ann, ", No: ")[end])
            catch
                missing
            end
        end
        #if ismissing(No[i])
            Annotations[i] = par.sub[i].ann
        #else
        #    Annotations[i] = join(split(par.sub[i].ann, ", ")[1:end-1], ", ")
        #end
    end  
    df = sort!(DataFrame(No = No, Name = Name, CAS = CAS, tR = tR, τR = τR, TR=TR, σR = σR, uR = uR, kR = kR, Annotations = Annotations, ), [:tR])
    Threads.@threads for i=1:n-1
        Res[i] = (df.tR[i+1] - df.tR[i])/(2*(df.τR[i+1] + df.τR[i]))
        Δs[i] = (df.tR[i+1] - df.tR[i])/(df.τR[i+1] - df.τR[i]) * log(df.τR[i+1]/df.τR[i])
    end
    df[!, :Res] = Res
    df[!, :Δs] = Δs 
    
    return select(df, [:No, :Name, :CAS, :tR, :τR, :TR, :σR, :uR, :kR, :Res, :Δs, :Annotations])
end

function peaklist_unthread(sol, par)
	n = length(par.sub)
    T = typeof(sol[1].u[end][1])
    # sol is solution from ODE system
    No = Array{Union{Missing, Int64}}(undef, n)
    Name = Array{String}(undef, n)
    CAS = Array{String}(undef, n)
    tR = Array{T}(undef, n)
    TR = Array{T}(undef, n)
    σR = Array{T}(undef, n)
    uR = Array{T}(undef, n)
    τR = Array{T}(undef, n)
    kR = Array{T}(undef, n)
    if T == Measurements.Measurement{Float64}
        Res = fill(NaN ± 0.0, n)
        Δs = fill(NaN ± 0.0, n)
    else
        Res = fill(NaN, n)
        Δs = fill(NaN, n)
    end
    Annotations = Array{String}(undef, n)
    for i=1:n
        Name[i] = par.sub[i].name
        CAS[i] = par.sub[i].CAS
        if sol[i].t[end]==par.col.L
            tR[i] = sol[i].u[end][1]
            TR[i] = par.prog.T_itp(par.col.L, tR[i]) - 273.15 
            uR[i] = 1/residency(par.col.L, tR[i], par.prog.T_itp, par.prog.Fpin_itp, par.prog.pout_itp, par.col.L, par.col.d, par.col.df, par.col.gas, par.sub[i].Tchar, par.sub[i].θchar, par.sub[i].ΔCp, par.sub[i].φ₀; ng=par.opt.ng, vis=par.opt.vis, control=par.opt.control, k_th=par.opt.k_th)
            τR[i] = sqrt(sol[i].u[end][2])
            σR[i] = τR[i]*uR[i]
            kR[i] = retention_factor(par.col.L, tR[i], par.prog.T_itp, par.col.d, par.col.df, par.sub[i].Tchar, par.sub[i].θchar, par.sub[i].ΔCp, par.sub[i].φ₀; k_th=par.opt.k_th)
        else
            tR[i] = NaN
            TR[i] = NaN
            uR[i] = NaN
            τR[i] = NaN
            σR[i] = NaN
            kR[i] = NaN
        end
        No[i] = try
            parse(Int,split(par.sub[i].ann, ", ")[end])
        catch
            try 
                parse(Int,split(par.sub[i].ann, ", No: ")[end])
            catch
                missing
            end
        end
        #if ismissing(No[i])
            Annotations[i] = par.sub[i].ann
        #else
        #    Annotations[i] = join(split(par.sub[i].ann, ", ")[1:end-1], ", ")
        #end
    end  
    df = sort!(DataFrame(No = No, Name = Name, CAS = CAS, tR = tR, τR = τR, TR=TR, σR = σR, uR = uR, kR = kR, Annotations = Annotations, ), [:tR])
    for i=1:n-1
        Res[i] = (df.tR[i+1] - df.tR[i])/(2*(df.τR[i+1] + df.τR[i]))
        Δs[i] = (df.tR[i+1] - df.tR[i])/(df.τR[i+1] - df.τR[i]) * log(df.τR[i+1]/df.τR[i])
    end
    df[!, :Res] = Res
    df[!, :Δs] = Δs 
    
    return select(df, [:No, :Name, :CAS, :tR, :τR, :TR, :σR, :uR, :kR, :Res, :Δs, :Annotations])
end


function peaklist_thread(sol, L, d, df, gas, T_itp, Fpin_itp, pout_itp, Name, CAS, ann, Tchar, θchar, ΔCp, φ₀, opt)
	n = length(Name)
    T = typeof(sol[1].u[end][1])
    # sol is solution from ODE system
    No = Array{Union{Missing, Int64}}(undef, n)
    tR = Array{T}(undef, n)
    TR = Array{T}(undef, n)
    σR = Array{T}(undef, n)
    uR = Array{T}(undef, n)
    τR = Array{T}(undef, n)
    kR = Array{T}(undef, n)
    if T == Measurements.Measurement{Float64}
        Res = fill(NaN ± 0.0, n)
        Δs = fill(NaN ± 0.0, n)
    else
        Res = fill(NaN, n)
        Δs = fill(NaN, n)
    end
    Threads.@threads for i=1:n
        if sol[i].t[end]==L
            tR[i] = sol[i].u[end][1]
            TR[i] = T_itp(L, tR[i]) - 273.15 
            uR[i] = 1/residency(L, tR[i], T_itp, Fpin_itp, pout_itp, L, d, df, gas, Tchar[i], θchar[i], ΔCp[i], φ₀[i]; ng=opt.ng, vis=opt.vis, control=opt.control, k_th=opt.k_th)
            τR[i] = sqrt(sol[i].u[end][2])
            σR[i] = τR[i]*uR[i]
            kR[i] = retention_factor(L, tR[i], T_itp, d, df, Tchar[i], θchar[i], ΔCp[i], φ₀[i]; k_th=opt.k_th)
        else
            tR[i] = NaN
            TR[i] = NaN
            uR[i] = NaN
            τR[i] = NaN
            σR[i] = NaN
            kR[i] = NaN
        end
        No[i] = try
            parse(Int,split(ann, ", ")[end])
        catch
            try 
                parse(Int,split(ann, ", No: ")[end])
            catch
                missing
            end
        end
    end  
    df = sort!(DataFrame(No = No, Name = Name, CAS = CAS, tR = tR, τR = τR, TR=TR, σR = σR, uR = uR, kR = kR, Annotations = ann, ), [:tR])
    Threads.@threads for i=1:n-1
        Res[i] = (df.tR[i+1] - df.tR[i])/(2*(df.τR[i+1] + df.τR[i]))
        Δs[i] = (df.tR[i+1] - df.tR[i])/(df.τR[i+1] - df.τR[i]) * log(df.τR[i+1]/df.τR[i])
    end
    df[!, :Res] = Res
    df[!, :Δs] = Δs 
    
    return select(df, [:No, :Name, :CAS, :tR, :τR, :TR, :σR, :uR, :kR, :Res, :Δs, :Annotations])
end

function peaklist_unthread(sol, L, d, df, gas, T_itp, Fpin_itp, pout_itp, Name, CAS, ann, Tchar, θchar, ΔCp, φ₀, opt)
	n = length(Name)
    T = typeof(sol[1].u[end][1])
    # sol is solution from ODE system
    No = Array{Union{Missing, Int64}}(undef, n)
    tR = Array{T}(undef, n)
    TR = Array{T}(undef, n)
    σR = Array{T}(undef, n)
    uR = Array{T}(undef, n)
    τR = Array{T}(undef, n)
    kR = Array{T}(undef, n)
    if T == Measurements.Measurement{Float64}
        Res = fill(NaN ± 0.0, n)
        Δs = fill(NaN ± 0.0, n)
    else
        Res = fill(NaN, n)
        Δs = fill(NaN, n)
    end
    for i=1:n
        if sol[i].t[end]==L
            tR[i] = sol[i].u[end][1]
            TR[i] = T_itp(L, tR[i]) - 273.15 
            uR[i] = 1/residency(L, tR[i], T_itp, Fpin_itp, pout_itp, L, d, df, gas, Tchar[i], θchar[i], ΔCp[i], φ₀[i]; ng=opt.ng, vis=opt.vis, control=opt.control, k_th=opt.k_th)
            τR[i] = sqrt(sol[i].u[end][2])
            σR[i] = τR[i]*uR[i]
            kR[i] = retention_factor(L, tR[i], T_itp, d, df, Tchar[i], θchar[i], ΔCp[i], φ₀[i]; k_th=opt.k_th)
        else
            tR[i] = NaN
            TR[i] = NaN
            uR[i] = NaN
            τR[i] = NaN
            σR[i] = NaN
            kR[i] = NaN
        end
        No[i] = try
            parse(Int,split(ann, ", ")[end])
        catch
            try 
                parse(Int,split(ann, ", No: ")[end])
            catch
                missing
            end
        end
    end  
    df = sort!(DataFrame(No = No, Name = Name, CAS = CAS, tR = tR, τR = τR, TR=TR, σR = σR, uR = uR, kR = kR, Annotations = ann, ), [:tR])
    for i=1:n-1
        Res[i] = (df.tR[i+1] - df.tR[i])/(2*(df.τR[i+1] + df.τR[i]))
        Δs[i] = (df.tR[i+1] - df.tR[i])/(df.τR[i+1] - df.τR[i]) * log(df.τR[i+1]/df.τR[i])
    end
    df[!, :Res] = Res
    df[!, :Δs] = Δs 
    
    return select(df, [:No, :Name, :CAS, :tR, :τR, :TR, :σR, :uR, :kR, :Res, :Δs, :Annotations])
end

"""
    peaklist(sol, peak, par)

Construct a DataFrame with the peak list of the solution `sol` and `peak` of
the simulation of the GC system defined by `par`. 

# Output

The peaklist DataFrame consists of the entrys:
* `No`: Number of the solute in the database.  
* `Name`: Name of the solute.
* `tR`: Retention time of the solute (in s).
* `τR`: Peak width of the solute (in s). 
* `TR`: Temperature of the end of the column at the retention time (in °C).
* `σR`: Band width of the solute at retention time (in m).
* `uR`: Solute velocity at retention time (in m/s).
* `kR`: Retention factor of the solute at retention time.
* `Res`: Resolution (4τ) between neighboring peaks.
* `Δs`: separation metric between neighboring peaks, assuming linear development of peak width `τR` between the peaks.
* `Annotations`: additional anotations, e.g. Source, categories, if available

# Examples

```julia
julia> pl = peaklist(sol, peak, par)
...
```    
"""
function peaklist(sol, peak, par)
	n = length(par.sub)
    T = typeof(sol[1].u[end])
    No = Array{Union{Missing, Int64}}(undef, n)
    Name = Array{String}(undef, n)
    CAS = Array{String}(undef, n)
    tR = Array{T}(undef, n)
    TR = Array{T}(undef, n)
    σR = Array{T}(undef, n)
    uR = Array{T}(undef, n)
    τR = Array{T}(undef, n)
    kR = Array{T}(undef, n)
    if T == Measurements.Measurement{Float64}
        Res = fill(NaN ± 0.0, n)
        Δs = fill(NaN ± 0.0, n)
    else
        Res = fill(NaN, n)
        Δs = fill(NaN, n)
    end
    Annotations = Array{String}(undef, n)
    #Threads.@threads for i=1:n
    for i=1:n
        Name[i] = par.sub[i].name
        CAS[i] = par.sub[i].CAS
        if sol[i].t[end]==par.col.L
            tR[i] = sol[i].u[end]
            TR[i] = par.prog.T_itp(par.col.L, tR[i]) - 273.15 
            uR[i] = 1/residency(par.col.L, tR[i], par.prog.T_itp, par.prog.Fpin_itp, par.prog.pout_itp, par.col.L, par.col.d, par.col.df, par.col.gas, par.sub[i].Tchar, par.sub[i].θchar, par.sub[i].ΔCp, par.sub[i].φ₀; ng=par.opt.ng, vis=par.opt.vis, control=par.opt.control, k_th=par.opt.k_th)
            τR[i] = sqrt(peak[i].u[end])
            σR[i] = τR[i]*uR[i]
            kR[i] = retention_factor(par.col.L, tR[i], par.prog.T_itp, par.col.d, par.col.df, par.sub[i].Tchar, par.sub[i].θchar, par.sub[i].ΔCp, par.sub[i].φ₀; k_th=par.opt.k_th)
        else
            tR[i] = NaN
            TR[i] = NaN
            uR[i] = NaN
            τR[i] = NaN
            σR[i] = NaN
            kR[i] = NaN
        end
        No[i] = try
            parse(Int,split(par.sub[i].ann, ", ")[end])
        catch
            missing
        end
        if ismissing(No[i])
            Annotations[i] = par.sub[i].ann
        else
            Annotations[i] = join(split(par.sub[i].ann, ", ")[1:end-1], ", ")
        end
    end  
    df = sort!(DataFrame(No = No, Name = Name, CAS = CAS, tR = tR, τR = τR, TR=TR, σR = σR, uR = uR, kR = kR, Annotations = Annotations, ), [:tR])
    #Threads.@threads for i=1:n-1
    for i=1:n-1
        Res[i] = (df.tR[i+1] - df.tR[i])/(2*(df.τR[i+1] + df.τR[i]))
        Δs[i] = (df.tR[i+1] - df.tR[i])/(df.τR[i+1] - df.τR[i]) * log(df.τR[i+1]/df.τR[i])
    end
    df[!, :Res] = Res 
    df[!, :Δs] = Δs   
    
    return select(df, [:No, :Name, :CAS, :tR, :τR, :TR, :σR, :uR, :kR, :Res, :Δs, :Annotations])
end

function peaklist(sol, peak, L, d, df, gas, T_itp, Fpin_itp, pout_itp, Name, CAS, ann, Tchar, θchar, ΔCp, φ₀, opt)
	n = length(Name)
    T = typeof(sol[1].u[end])
    No = Array{Union{Missing, Int64}}(undef, n)
    tR = Array{T}(undef, n)
    TR = Array{T}(undef, n)
    σR = Array{T}(undef, n)
    uR = Array{T}(undef, n)
    τR = Array{T}(undef, n)
    kR = Array{T}(undef, n)
    if T == Measurements.Measurement{Float64}
        Res = fill(NaN ± 0.0, n)
        Δs = fill(NaN ± 0.0, n)
    else
        Res = fill(NaN, n)
        Δs = fill(NaN, n)
    end
    #Threads.@threads for i=1:n
    for i=1:n
        if sol[i].t[end]==L
            tR[i] = sol[i].u[end]
            TR[i] = T_itp(L, tR[i]) - 273.15 
            uR[i] = 1/residency(L, tR[i], T_itp, Fpin_itp, pout_itp, L, d, df, gas, Tchar[i], θchar[i], ΔCp[i], φ₀[i]; ng=opt.ng, vis=opt.vis, control=opt.control, k_th=opt.k_th)
            τR[i] = sqrt(peak[i].u[end])
            σR[i] = τR[i]*uR[i]
            kR[i] = retention_factor(L, tR[i], T_itp, d, df, Tchar[i], θchar[i], ΔCp[i], φ₀[i]; k_th=opt.k_th)
        else
            tR[i] = NaN
            TR[i] = NaN
            uR[i] = NaN
            τR[i] = NaN
            σR[i] = NaN
            kR[i] = NaN
        end
        No[i] = try
            parse(Int,split(ann, ", ")[end])
        catch
            missing
        end
    end  
    df = sort!(DataFrame(No = No, Name = Name, CAS = CAS, tR = tR, τR = τR, TR=TR, σR = σR, uR = uR, kR = kR, Annotations = ann, ), [:tR])
    #Threads.@threads for i=1:n-1
    for i=1:n-1
        Res[i] = (df.tR[i+1] - df.tR[i])/(2*(df.τR[i+1] + df.τR[i]))
        Δs[i] = (df.tR[i+1] - df.tR[i])/(df.τR[i+1] - df.τR[i]) * log(df.τR[i+1]/df.τR[i])
    end
    df[!, :Res] = Res 
    df[!, :Δs] = Δs   
    
    return select(df, [:No, :Name, :CAS, :tR, :τR, :TR, :σR, :uR, :kR, :Res, :Δs, :Annotations])
end


"""
    sol_extraction(sol, par)

Extract the points z=t, t=u1, τ²=u2 from the solution `sol` of
the ODE system of the GC system defined by `par` and exports them in a DataFrame.

# Examples

```julia
df_sol = sol_extraction(sol, par)
...
```    
"""
function sol_extraction(sol, par)
    # extract the points z=t, t=u1, τ²=u2 from the solution of
    # the ODE system
	n = length(par.sub)
    sol_z = Array{typeof(sol[1].t)}(undef, n)
    sol_t = Array{typeof(sol[1].u[1])}(undef, n)
    sol_τ² = Array{typeof(sol[1].u[1])}(undef, n)
    solutes = Array{String}(undef, n)
    for i=1:n
        sol_z[i] = sol[i].t
        temp_t = Array{typeof(sol[1].u[end][1])}(undef, length(sol[i].t))
        temp_τ² = Array{typeof(sol[1].u[end][2])}(undef, length(sol[i].t))
        for j=1:length(sol[i].t)
                temp_t[j] = sol[i].u[j][1]
                temp_τ²[j] = sol[i].u[j][2]
        end
        sol_t[i] = temp_t
        sol_τ²[i] = temp_τ²
        solutes[i] = par.sub[i].name
    end
    df_sol = DataFrame(name=solutes, z=sol_z, t=sol_t, τ²=sol_τ²)
    return df_sol
end

"""
    sol_extraction(sol, peak, par)

Extract the points z_t=sol.t, t=sol.u, z_τ²=peak.t and τ²=peak.u from the
solution `sol` and `peak` of the ODEs of the GC system defined by `par` and exports them in a DataFrame.

# Examples

```julia
df_sol = sol_extraction(sol, peak, par)
...
```    
"""
function sol_extraction(sol, peak, par)
    # extract the points z=t, t=u, from the solution of
    # the first ODE (sol_tz)
    # and the points z=t, τ²=u, from the solution of 
    # the second ODE (peak_τz)
	n = length(par.sub)
    sol_z = Array{typeof(sol[1].t)}(undef, n)
    sol_t = Array{typeof(sol[1].u)}(undef, n)
    peak_z = Array{typeof(peak[1].t)}(undef, n)
    peak_τ² = Array{typeof(peak[1].u)}(undef, n)
    solutes = Array{String}(undef, n)
    for i=1:n
        sol_z[i] = sol[i].t
        sol_t[i] = sol[i].u
        peak_z[i] = peak[i].t
        peak_τ²[i] = peak[i].u
		solutes[i] = par.sub[i].name
	end
    df_sol = DataFrame(name=solutes, z_t=sol_z, t=sol_t, z_τ²=peak_z, τ²=peak_τ²)
    return df_sol
end


end # module
