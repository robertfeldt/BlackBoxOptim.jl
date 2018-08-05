"""
Embedding operator that randomly samples
between parent's value and the nearest parameter boundary
to get the new valid value if target's parameter is out-of-bounds.
"""
struct RandomBound{S<:SearchSpace} <: EmbeddingOperator
    searchSpace::S

    RandomBound(searchSpace::S) where {S<:SearchSpace} = new{S}(searchSpace)
end

# outer ctors
RandomBound(dimBounds::Vector{ParamBounds}) = RandomBound(RangePerDimSearchSpace(dimBounds))

search_space(rb::RandomBound) = rb.searchSpace

@inline function random_bound(t::Real, r::Real, l::Real, u::Real)
    if t < l
        t = l + rand() * (r - l)
    elseif t > u
        t = u + rand() * (r - u)
    end
    @assert l <= t <= u "target=$(t) is out of [$l, $u]" # if r is out of [l, u]
    return t
end

function apply!(eo::RandomBound, target::AbstractIndividual, ref::AbstractIndividual)
    length(target) == length(ref) == numdims(eo.searchSpace) ||
        throw(ArgumentError("Dimensions of problem/individuals do not match"))
    ss = search_space(eo)
    @inbounds target .= random_bound.(target, ref, mins(ss), maxs(ss))
    return target
end

apply!(eo::RandomBound, target::AbstractIndividual, pop, refIndex::Int) =
    apply!(eo, target, viewer(pop, refIndex))

function apply!(eo::RandomBound, target::AbstractIndividual,
                pop, parentIndices::AbstractVector{Int})
    @assert length(parentIndices) == 1
    apply!(eo, target, pop, parentIndices[1])
end
