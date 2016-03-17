"""
  Embedding operator that randomly samples
  between parent's value and the nearest parameter boundary
  to get the new valid value if target's parameter is out-of-bounds.
"""
type RandomBound{S<:SearchSpace} <: EmbeddingOperator
    searchSpace::SearchSpace
end

# outer ctors
RandomBound{S<:SearchSpace}(searchSpace::S) = RandomBound{S}(searchSpace)
RandomBound(dimBounds::Vector{ParamBounds}) = RandomBound(RangePerDimSearchSpace(dimBounds))

search_space(rb::RandomBound) = rb.searchSpace

function apply!(eo::RandomBound, target::Individual, ref::AbstractVector)
  length(target) == length(ref) == numdims(eo.searchSpace) || throw(ArgumentError("Dimensions of problem/individuals do not match"))
  ssmins = mins(eo.searchSpace)
  ssmaxs = maxs(eo.searchSpace)

  for i in 1:length(target)
    l, u = ssmins[i], ssmaxs[i]

    if target[i] < l
      target[i] = l + rand() * (ref[i] - l)
    elseif target[i] > u
      target[i] = u + rand() * (ref[i] - u)
    end
    @assert l <= target[i] <= u "target[$i]=$(target[i]) is out of [$l, $u]"
  end
  return target
end

apply!(eo::RandomBound, target::Individual, pop, refIndex::Int) =
  apply!(eo, target, view(pop, refIndex))

function apply!(eo::RandomBound, target::Individual, pop, parentIndices::Vector{Int})
  @assert length(parentIndices) == 1
  apply!(eo, target, pop, parentIndices[1])
end
