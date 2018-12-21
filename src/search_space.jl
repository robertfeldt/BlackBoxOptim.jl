"""
A `SearchSpace` defines the valid candidate points that could be
considered in a search/optimization.
The base abstract class has very few restrictions
and can allow varying number of dimensions etc.
"""
abstract type SearchSpace end

"""
`SearchSpace` with a fixed finite number of dimensions.
Applicable to the vast majority of problems.
"""
abstract type FixedDimensionSearchSpace <: SearchSpace end

"""
Fixed-dimensional space, each dimension has a continuous range of valid values.
"""
abstract type ContinuousSearchSpace <: FixedDimensionSearchSpace end

"""
The point of the `SearchSpace`.

The abstract type. It allows different actual implementations to be used,
e.g `Vector` or `SubArray`.
"""
const AbstractIndividual = AbstractVector{Float64}

"""
The point of the `SearchSpace`.

The concrete type that could be used for storage.
"""
const Individual = Vector{Float64}

"""
The valid range of values for a specific dimension in a `SearchSpace`.
"""
const ParamBounds = Tuple{Float64,Float64}

"""
Get the range of valid values for a specific dimension.
"""
range_for_dim(ss::ContinuousSearchSpace, i) = (mins(ss)[i], maxs(ss)[i])

ranges(ss::ContinuousSearchSpace) = collect(zip(mins(ss), maxs(ss)))

"""
Generate `numIndividuals` individuals by random sampling in the search space.
"""
rand_individuals(ss::ContinuousSearchSpace, numIndividuals) =
    mins(ss) .+ deltas(ss) .* rand(numdims(ss), numIndividuals)

"""
Generate `numIndividuals` individuals by latin hypercube sampling (LHS).
This should be the default way to create the initial population.
"""
rand_individuals_lhs(ss::ContinuousSearchSpace, numIndividuals) =
    Utils.latin_hypercube_sampling(mins(ss), maxs(ss), numIndividuals)

"""
Generate one random candidate.
"""
rand_individual(ss::ContinuousSearchSpace) =
    mins(ss) .+ deltas(ss) .* rand(numdims(ss))

"""
Check if given individual lies in the given search space.
"""
function Base.in(ind::AbstractIndividual, ss::ContinuousSearchSpace)
    @assert length(ind) == numdims(ss)
    @inbounds for i in eachindex(ind)
        (mins(ss)[i] <= ind[i] <= maxs(ss)[i]) || return false
    end
    return true
end

"""
`SearchSpace` defined by a range of valid values for each dimension.
"""
struct RangePerDimSearchSpace <: ContinuousSearchSpace
    # We save the ranges as individual mins, maxs and deltas for faster access later.
    mins::Vector{Float64}
    maxs::Vector{Float64}
    deltas::Vector{Float64}

    function RangePerDimSearchSpace(ranges)
        mins = getindex.(ranges, 1)
        maxs = getindex.(ranges, 2)
        new(mins, maxs, maxs .- mins)
    end

    RangePerDimSearchSpace(mins, maxs) = new(mins, maxs, maxs .- mins)
end

mins(ss::RangePerDimSearchSpace) = ss.mins
maxs(ss::RangePerDimSearchSpace) = ss.maxs
deltas(ss::RangePerDimSearchSpace) = ss.deltas
numdims(ss::RangePerDimSearchSpace) = length(mins(ss))

diameters(ss::RangePerDimSearchSpace) = deltas(ss)

Base.:(==)(a::RangePerDimSearchSpace, b::RangePerDimSearchSpace) =
    numdims(a) == numdims(b) && (a.mins == b.mins) && (a.maxs == b.maxs)

"""
Create `RangePerDimSearchSpace` with given number of dimensions
and given range of valid values for each dimension.
"""
symmetric_search_space(numdims, range=(0.0, 1.0)) = RangePerDimSearchSpace(fill(range, numdims))

"""
Projects a given point onto the search space coordinate-wise.
"""
feasible(v::AbstractIndividual, ss::RangePerDimSearchSpace) = clamp.(v, mins(ss), maxs(ss))

# concatenates two range-based search spaces
Base.vcat(ss1::RangePerDimSearchSpace, ss2::RangePerDimSearchSpace) =
    RangePerDimSearchSpace(vcat(mins(ss1), mins(ss2)),
                           vcat(maxs(ss1), maxs(ss2)))

"""
0-dimensional search space.
Could be used as a placeholder for optional `SearchSpace` parameters.
"""
const ZERO_SEARCH_SPACE = RangePerDimSearchSpace(Vector{Float64}(), Vector{Float64}())
