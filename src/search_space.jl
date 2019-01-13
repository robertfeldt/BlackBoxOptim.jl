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
I.e. `dimmin(ss, i)` ≤ `x[i]` ≤ `dimmax(ss, i)` for each dimension `i`.
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
    dimmin(ss::SearchSpace, [i::Integer])

A minimal valid value for `i`-th dimension of `ss`, or a vector of
minimal valid values for each dimension of `ss` if no `i` was given.
"""
dimmin(ss::RectSearchSpace) = ss.dimmin

dimmin(ss::RectSearchSpace, i::Integer) = ss.dimmin[i]

"""
    dimmax(ss::SearchSpace, [i::Integer])

A maximal valid value for `i`-th dimension of `ss`, or a vector of
maximal valid values for each dimension of `ss` if no `i` was given.
"""
dimmax(ss::RectSearchSpace) = ss.dimmax

dimmax(ss::RectSearchSpace, i::Integer) = ss.dimmax[i]

"""
    dimdelta(ss::SearchSpace, [i::Integer])

A delta of maximal and minimal valid values for `i`-th dimension of `ss`
(*diameter* of `i`-th dimension), or a vector of deltas for each dimension
if no `i` given.
"""
dimdelta(ss::RectSearchSpace) = ss.dimdelta

dimdelta(ss::RectSearchSpace, i::Integer) = ss.dimdelta[i]

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
dimrange(ss::RectSearchSpace, i::Integer) = (dimmin(ss, i), dimmax(ss, i))
@deprecate range_for_dim(ss, i) dimrange(ss, i)

dimrange(ss::RectSearchSpace) = tuple.(dimmin(ss), dimmax(ss))
@deprecate ranges(ss) dimrange(ss)

"""
    in(ind::AbstractIndividual, ss::SearchSpace)

Check if given individual lies in the given search space.
"""
function Base.in(ind::AbstractIndividual, ss::RectSearchSpace)
    # FIXME checks for the discrete space?
    @assert length(ind) == numdims(ss)
    @inbounds for i in eachindex(ind)
        (dimmin(ss, i) <= ind[i] <= dimmax(ss, i)) || return false
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
Projects a given point onto the search space coordinate-wise.
"""
feasible(v::AbstractIndividual, ss::RectSearchSpace) =
    clamp.(v, dimmin(ss), dimmax(ss))

# concatenates two range-based search spaces
Base.vcat(ss1::ContinuousRectSearchSpace,
          ss2::ContinuousRectSearchSpace) =
    ContinuousRectSearchSpace(vcat(dimmin(ss1), dimmin(ss2)),
                              vcat(dimmax(ss1), dimmax(ss2)))

# all dimensions are continuous
dimdigits(ss::ContinuousRectSearchSpace, i::Integer) = -1
dimdigits(ss::ContinuousRectSearchSpace) = fill(-1, numdims(ss))

"""
0-dimensional search space.
Could be used as a placeholder for optional `SearchSpace` parameters.
"""
const ZERO_SEARCH_SPACE = ContinuousRectSearchSpace(Vector{Float64}(), Vector{Float64}())

"""
`RectSearchSpace` that allows both continuous and discrete dimensions.
Discrete dimensions are defined by the number of digits (`dimdigits(ss)`) after decimal point (i.e. precision)

If `dimdigits(ss, i)` is negative, `i`-th dimension is considered continuous.
"""
struct MixedPrecisionRectSearchSpace <: RectSearchSpace
    dimmin::Vector{Float64}
    dimmax::Vector{Float64}
    dimdelta::Vector{Float64}
    dimdigits::Vector{Int}

    function MixedPrecisionRectSearchSpace(dimmin::AbstractVector,
                                           dimmax::AbstractVector,
                                           dimdigits::AbstractVector{<:Integer})
        length(dimmin) == length(dimmax) == length(dimdigits) ||
            throw(DimensionMismatch("dimmin, dimmax and dimdigits should have the same length"))
        dmin = Float64[dimdigits[i] >= 0 ? round(dimmin[i], digits=dimdigits[i]) : dimmin[i] for i in eachindex(dimdigits)]
        dmax = Float64[dimdigits[i] >= 0 ? round(dimmax[i], digits=dimdigits[i]) : dimmax[i] for i in eachindex(dimdigits)]
        all(xy -> xy[1] <= xy[2], zip(dmin, dmax)) ||
            throw(ArgumentError("discretized dimmin should not exceed dimmax"))
        new(dmin, dmax, dmax .- dmin,
            copyto!(Vector{Int}(undef, length(dimdigits)), dimdigits))
    end
end

MixedPrecisionRectSearchSpace(dimranges::AbstractVector{<:Tuple{<:Number, <:Number}},
                              dimdigits::AbstractVector{<:Integer}) =
    MixedPrecisionRectSearchSpace(getindex.(dimranges, 1), getindex.(dimranges, 2), dimdigits)

dimdigits(ss::MixedPrecisionRectSearchSpace, i::Integer) = ss.dimdigits[i]
dimdigits(ss::MixedPrecisionRectSearchSpace) = ss.dimdigits

feasible(x::Real, u::Real, v::Real, dimdigits::Integer) =
    clamp(dimdigits >= 0 ? round(x, digits=dimdigits) : x, u, v)

"""
Projects a given point onto the search space coordinate-wise.
"""
feasible(v::AbstractIndividual, ss::MixedPrecisionRectSearchSpace) =
    feasible.(v, dimmin(ss), dimmax(ss), dimdigits(ss))

# concatenates two rectangular search spaces
# one of the should be MixedPrecisionRectSearchSpace, otherwise
# ContinuousRectSearchSpace-only method would be called
Base.vcat(ss1::RectSearchSpace, ss2::RectSearchSpace) =
    MixedPrecisionRectSearchSpace(vcat(dimmin(ss1), dimmin(ss2)),
                                  vcat(dimmax(ss1), dimmax(ss2)),
                                  vcat(dimdigits(ss1), dimdigits(ss2)))

"""
    RectSearchSpace(dimranges::AbstractVector; [dimdigits = nothing])

Create `RectSearchSpace` with given range of valid values for each dimension and,
if specified, `dimdigits` precision for each dimension.
Returns `MixedPrecisionRectSearchSpace` if there's at least one dimension
with specified precision (dimdigits[i] ≥ 0), otherwise `ContinuousRectSearchSpace`.
"""
RectSearchSpace(dimranges::AbstractVector;
                dimdigits::Union{AbstractVector{<:Integer}, Nothing} = nothing) =
    dimdigits === nothing || all(ddigits -> ddigits < 0, dimdigits) ?
        ContinuousRectSearchSpace(dimranges) :
        MixedPrecisionRectSearchSpace(dimranges, dimdigits)

"""
Create `RectSearchSpace` with given number of dimensions
and given range of valid values for each dimension.
"""
RectSearchSpace(numdims::Integer, range=(0.0, 1.0); dimdigits::Union{Integer, Nothing} = nothing) =
    RectSearchSpace(fill(range, numdims);
                    dimdigits = dimdigits !== nothing && dimdigits >= 0 ? fill(dimdigits, numdims) : nothing)

@deprecate symmetric_search_space(numdims, range=(0.0, 1.0); dimdigits=nothing) RectSearchSpace(numdims, range, dimdigits=dimdigits)

# round x to fit ss requirements
# by default it does nothing
_round!(x::Union{AbstractVector, AbstractMatrix}, ss::RectSearchSpace) = x

function _round!(x::AbstractVector, ss::MixedPrecisionRectSearchSpace)
    length(x) == numdims(ss) ||
        throw(DimensionMismatch("Vector length should match search space dimensions"))
    ddigits = dimdigits(ss)
    @inbounds for j in eachindex(x)
        if ddigits[j] >= 0
            x[j] = round(x[j], digits=ddigits[j])
        end
    end
    return x
end

function _round!(x::AbstractMatrix, ss::MixedPrecisionRectSearchSpace)
    size(x, 1) == numdims(ss) ||
        throw(DimensionMismatch("Matrix rows should match search space dimensions"))
    ddigits = dimdigits(ss)
    @inbounds for i in axes(x, 2)
        for j in axes(x, 1)
            if ddigits[j] >= 0
                x[j, i] = round(x[j, i], digits=ddigits[j])
            end
        end
    end
    return x
end

"""
Generate one random candidate.
"""
rand_individual(ss::RectSearchSpace) =
    _round!(rand(numdims(ss)) .* dimdelta(ss) .+ dimmin(ss), ss)

"""
    rand_individuals(ss, n; [method=:latin_hypercube])

Generate `n` individuals by randomly sampling in the `ss` search space using
specified sampling `method`. The supported methods are:
 * `uniform`: uniform independent sampling for each dimension
 * `latin_hypercube` (the default): use latin hypercube sampling method
"""
function rand_individuals(ss::RectSearchSpace, n::Integer; method::Symbol=:latin_hypercube)
    if method == :uniform
        return _round!(rand(numdims(ss), n) .* dimdelta(ss) .+ dimmin(ss), ss)
    elseif method == :latin_hypercube
        return _round!(Utils.latin_hypercube_sampling(dimmin(ss), dimmax(ss), n), ss)
    else
        throw(ArgumentError("Unknown sampling method \"$method\""))
    end
end

@deprecate rand_individuals_lhs(ss::RectSearchSpace, n::Integer) rand_individuals(ss, n, method=:latin_hypercube)
