# Embedding operator operator that randomly samples between parent's value and the parameter boundary
# to get the new valid value if target's parameter is out-of-bounds
type RandomBound{S<:SearchSpace} <: EmbeddingOperator
    searchSpace::SearchSpace
end

# outer ctors
RandomBound{S<:SearchSpace}(searchSpace::S) = RandomBound{S}(searchSpace)
RandomBound(dimBounds::Vector{ParamBounds}) = RandomBound(RangePerDimSearchSpace(dimBounds))

function apply!(eo::RandomBound, target::Individual, pop, parentIndices)
  @assert length(parentIndices) == 1
  ssmins = mins(eo.searchSpace)
  ssmaxs = maxs(eo.searchSpace)

  parentIx = parentIndices[1]
  for i in 1:length(target)
    l, u = ssmins[i], ssmaxs[i]

    if target[i] < l
      target[i] = l + rand() * (pop[i, parentIx] - l)
    elseif target[i] > u
      target[i] = u + rand() * (pop[i, parentIx] - u)
    end
  end
  return target
end

# 1-individual population
apply!(eo::RandomBound, target::Individual, example::Individual) = apply!(eo, target, example, [1])
