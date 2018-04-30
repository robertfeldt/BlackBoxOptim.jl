"""
A base abstract type for `OptimizationProblem` search space.
A concrete `SearchSpace` subtype specifies the valid
candidate points that could be considered in a search/optimization.
"""
abstract type SearchSpace end

"""
A base abstract type for search spaces with a fixed finite number of dimensions.
Applicable to the vast majority of problems.
"""
abstract type FixedDimensionSearchSpace <: SearchSpace end

"""
A `SearchSpace` with `N`-dimensional rectangle as a set of valid points.
I.e. `mins(ss)[i]` ≤ `x[i]` ≤ `maxs(ss)[i]` for each dimension `i`.
"""
abstract type RectSearchSpace <: FixedDimensionSearchSpace end

Base.@deprecate_binding ContinuousSearchSpace RectSearchSpace

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
The valid range of values for a specific dimension in a `RectSearchSpace`.
"""
const ParamBounds = Tuple{Float64,Float64}

numdims(ss::RectSearchSpace) = length(mins(ss))

mins(ss::RectSearchSpace) = ss.mins
maxs(ss::RectSearchSpace) = ss.maxs
deltas(ss::RectSearchSpace) = ss.deltas
diameters(ss::RectSearchSpace) = deltas(ss)

"""
Get the range of valid values for a specific dimension.
"""
range_for_dim(ss::RectSearchSpace, i::Integer) = (mins(ss)[i], maxs(ss)[i])
ranges(ss::RectSearchSpace) = collect(zip(dimmins(ss), dimmaxs(ss)))

"""
Check if given individual lies in the given search space.
"""
function Base.in(ind::AbstractIndividual, ss::RectSearchSpace)
    @assert length(ind) == numdims(ss)
    @inbounds for i in eachindex(ind)
        (mins(ss)[i] <= ind[i] <= maxs(ss)[i]) || return false
    end
    return true
end

"""
`SearchSpace` defined by a continuous range of valid values for each dimension.
"""
struct ContinuousRectSearchSpace <: RectSearchSpace
    # We save the ranges as individual mins, maxs and deltas for faster access later.
    mins::Vector{Float64}
    maxs::Vector{Float64}
    deltas::Vector{Float64}

    function ContinuousRectSearchSpace(
        mins::AbstractVector{<:Real},
        maxs::AbstractVector{<:Real}
    )
        length(mins) == length(maxs) ||
            throw(DimensionMismatch("mins and maxs should have the same length"))
        all(xy -> xy[1] <= xy[2], zip(mins, maxs)) ||
            throw(ArgumentError("mins should not exceed maxs"))
        new(mins, maxs, maxs .- mins)
    end
end

ContinuousRectSearchSpace(ranges) =
    ContinuousRectSearchSpace(getindex.(ranges, 1), getindex.(ranges, 2))

Base.:(==)(a::ContinuousRectSearchSpace,
           b::ContinuousRectSearchSpace) =
    numdims(a) == numdims(b) && (a.mins == b.mins) && (a.maxs == b.maxs)

"""
Generate one random candidate.
"""
rand_individual(ss::ContinuousRectSearchSpace) =
    mins(ss) .+ deltas(ss) .* rand(numdims(ss))

"""
Generate `n` individuals by random sampling in the search space.
"""
rand_individuals(ss::ContinuousRectSearchSpace, n::Integer) =
    mins(ss) .+ deltas(ss) .* rand(numdims(ss), n)

"""
Generate `n` individuals by latin hypercube sampling (LHS).
This should be the default way to create the initial population.
"""
rand_individuals_lhs(ss::ContinuousRectSearchSpace, n::Integer) =
    Utils.latin_hypercube_sampling(mins(ss), maxs(ss), n)

"""
Projects a given point onto the search space coordinate-wise.
"""
feasible(v::AbstractIndividual, ss::ContinuousRectSearchSpace) =
    clamp.(v, mins(ss), maxs(ss))

# concatenates two range-based search spaces
Base.vcat(ss1::ContinuousRectSearchSpace,
          ss2::ContinuousRectSearchSpace) =
    ContinuousRectSearchSpace(vcat(mins(ss1), mins(ss2)),
                              vcat(maxs(ss1), maxs(ss2)))

"""
0-dimensional search space.
Could be used as a placeholder for optional `SearchSpace` parameters.
"""
const ZERO_SEARCH_SPACE = ContinuousRectSearchSpace(Vector{Float64}(), Vector{Float64}())

"""
Create `RectSearchSpace` with given range of valid values for each dimension.
"""
RectSearchSpace(ranges::AbstractVector) =
    ContinuousRectSearchSpace(ranges)

"""
Create `RectSearchSpace` with given number of dimensions
and given range of valid values for each dimension.
"""
symmetric_search_space(numdims::Integer, range=(0.0, 1.0)) =
    RectSearchSpace(fill(range, numdims))
