"""
    Simplex Crossover (SPX).

    `ϵ>0` controls how the original simplex is inflated, `ϵ=1` means no
    inflation.

    See Tsutsui, Yamamura & Higuchi "Multi-parent recombination with simplex crossover in real coded genetic algorithms", 1999,
    Proc. of the Genetic and Evolutionary Computation Conference
"""
immutable SimplexCrossover{NP} <: CrossoverOperator{NP,1}
    ϵ::Float64      # inflation rate

    SimplexCrossover(ϵ::Number = Math.sqrt(NP + 1)) = new(ϵ)
    SimplexCrossover(params::Parameters) = new(params[:SPX_ϵ])
end

const SPX_DefaultOptions = ParamsDict(
  :SPX_ϵ => 1.1 # how much to inflate
)

#masscenter(pop, parentIndices) = mapslices(mean, pop[:, parentIndices], 2)

function apply!{NP}(xover::SimplexCrossover{NP}, target::Individual, targetIndex::Int, pop, parentIndices)
    @assert length(parentIndices) == NP

    offset = zeros(target)
    parentSum = copy(view(pop, parentIndices[1]))
    for i in 1:(NP-1)
        s = rand()^(1/i)
        parent = view(pop, parentIndices[i])
        parentNext = view(pop, parentIndices[i+1])
        @inbounds for j in 1:numdims(pop)
            offset[j] = s * (xover.ϵ*(parent[j] - parentNext[j]) + offset[j])
            parentSum[j] += parentNext[j]
        end
    end
    s = (1-xover.ϵ)/NP
    parentLast = view(pop, parentIndices[end])
    @inbounds for j in 1:numdims(pop)
        target[j] = s*parentSum[j] + xover.ϵ*parentLast[j] + offset[j]
    end
    return target
end
