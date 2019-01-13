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
I.e. `dimmin(ss)[i]` ≤ `x[i]` ≤ `dimmax(ss)[i]` for each dimension `i`.
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
The valid range of values for a specific dimension of `SearchSpace`.
"""
const ParamBounds = Tuple{Float64,Float64}

numdims(ss::RectSearchSpace) = length(dimmin(ss))

"""
    dimmin(ss::SearchSpace)

A vector of minimal valid values for each dimension of `ss`.
"""
dimmin(ss::RectSearchSpace) = ss.dimmin

"""
    dimmax(ss::SearchSpace)

A vector of maximal valid values for each dimension of `ss`.
"""
dimmax(ss::RectSearchSpace) = ss.dimmax

"""
    dimdelta(ss::SearchSpace)

A vector of deltas between maximal and minimal valid values
for each dimension of `ss`.
"""
dimdelta(ss::RectSearchSpace) = ss.dimdelta

@deprecate mins(ss) dimmin(ss)
@deprecate maxs(ss) dimmax(ss)
@deprecate deltas(ss) dimdelta(ss)
@deprecate diameters(ss) dimdelta(ss)

"""
    dimrange(ss::SearchSpace, [i::Integer])

Gets a `ParamsRange` tuple of minimal and maximal valid values for
`i`-th dimension of `ss`, or a vector of `ParamsRange` tuples
for each dimension if no `i` given.
"""
dimrange(ss::RectSearchSpace, i::Integer) = (dimmin(ss)[i], dimmax(ss)[i])
@deprecate range_for_dim(ss, i) dimrange(ss, i)

dimrange(ss::RectSearchSpace) = tuple.(dimmin(ss), dimmax(ss))
@deprecate ranges(ss) dimrange(ss)

"""
Check if given individual lies in the given search space.
"""
function Base.in(ind::AbstractIndividual, ss::RectSearchSpace)
    @assert length(ind) == numdims(ss)
    @inbounds for i in eachindex(ind)
        (dimmin(ss)[i] <= ind[i] <= dimmax(ss)[i]) || return false
    end
    return true
end

"""
`SearchSpace` defined by a continuous range of valid values for each dimension.
"""
struct ContinuousRectSearchSpace <: RectSearchSpace
    dimmin::Vector{Float64}    # minimal valid value per dimension
    dimmax::Vector{Float64}    # maximal valid value per dimension
    dimdelta::Vector{Float64}  # delta/diameter (dimmax-dimmin) per dimension

    function ContinuousRectSearchSpace(
        dimmin::AbstractVector{<:Real},
        dimmax::AbstractVector{<:Real}
    )
        length(dimmin) == length(dimmax) ||
            throw(DimensionMismatch("dimmin and dimmax should have the same length"))
        all(xy -> xy[1] <= xy[2], zip(dimmin, dimmax)) ||
            throw(ArgumentError("dimmin should not exceed dimmax"))
        new(dimmin, dimmax, dimmax .- dimmin)
    end
end

ContinuousRectSearchSpace(ranges) =
    ContinuousRectSearchSpace(getindex.(ranges, 1), getindex.(ranges, 2))

Base.:(==)(a::ContinuousRectSearchSpace,
           b::ContinuousRectSearchSpace) =
    (numdims(a) == numdims(b)) &&
    (dimmin(a) == dimmin(b)) &&
    (dimmax(a) == dimmax(b))

"""
Generate one random candidate.
"""
rand_individual(ss::ContinuousRectSearchSpace) =
    dimmin(ss) .+ dimdelta(ss) .* rand(numdims(ss))

"""
Generate `n` individuals by random sampling in the search space.
"""
rand_individuals(ss::ContinuousRectSearchSpace, n::Integer) =
    dimmin(ss) .+ dimdelta(ss) .* rand(numdims(ss), n)

"""
Generate `n` individuals by latin hypercube sampling (LHS).
This should be the default way to create the initial population.
"""
rand_individuals_lhs(ss::ContinuousRectSearchSpace, n::Integer) =
    Utils.latin_hypercube_sampling(dimmin(ss), dimmax(ss), n)

"""
Projects a given point onto the search space coordinate-wise.
"""
feasible(v::AbstractIndividual, ss::RectSearchSpace) =
    clamp.(v, dimmin(ss), dimmax(ss))

# concatenates two range-based search spaces
Base.vcat(ss1::ContinuousRectSearchSpace,
          ss2::ContinuousRectSearchSpace) =
    ContinuousRectSearchSpace(vcat(dimmin(ss1), dimmin(ss2)),
                              vcat(dimmax(ss1), dimmax(ss2)))

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
