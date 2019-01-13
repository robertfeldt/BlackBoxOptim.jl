"""
Embedding operator that randomly samples
between parent's value and the nearest parameter boundary
to get the new valid value if target's parameter is out-of-bounds.
"""
struct RandomBound{S<:RectSearchSpace} <: EmbeddingOperator
    search_space::S

    RandomBound(search_space::S) where {S<:RectSearchSpace} = new{S}(search_space)
end

search_space(rb::RandomBound) = rb.search_space

function apply!(eo::RandomBound, target::AbstractIndividual, ref::AbstractIndividual)
    length(target) == length(ref) == numdims(eo.search_space) ||
        throw(ArgumentError("Dimensions of problem/individuals do not match"))
    ss = search_space(eo)
    ssmins = dimmin(ss)
    ssmaxs = dimmax(ss)

    @inbounds for i in eachindex(target)
        l, u = ssmins[i], ssmaxs[i]

        if target[i] < l
            target[i] = l + rand() * (ref[i]-l)
        elseif target[i] > u
            target[i] = u + rand() * (ref[i]-u)
        else # continuous range doesn't need further checks
            (ss isa MixedPrecisionRectSearchSpace) || continue
        end
        if (ss isa MixedPrecisionRectSearchSpace) && (dimdigits(ss, i) >= 0)
            target[i] = round(target[i], digits=dimdigits(ss, i))
        end
        @assert l <= target[i] <= u "target[$i]=$(target[i]) is out of [$l, $u]"
    end
    return target
end

apply!(eo::RandomBound, target::AbstractIndividual, pop, refIndex::Int) =
    apply!(eo, target, viewer(pop, refIndex))

function apply!(eo::RandomBound, target::AbstractIndividual,
                pop, parentIndices::AbstractVector{Int})
    @assert length(parentIndices) == 1
    apply!(eo, target, pop, parentIndices[1])
end
