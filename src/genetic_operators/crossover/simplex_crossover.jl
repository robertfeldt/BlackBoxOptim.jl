"""
Simplex Crossover (SPX).

`ϵ>0` controls how the original simplex is inflated, `ϵ=1` means no
inflation.

See Tsutsui, Yamamura & Higuchi "Multi-parent recombination with simplex crossover in real coded genetic algorithms", 1999,
Proc. of the Genetic and Evolutionary Computation Conference
"""
struct SimplexCrossover{NP} <: CrossoverOperator{NP,1}
    ϵ::Float64      # inflation rate

    SimplexCrossover{NP}(ϵ::Number = Math.sqrt(NP + 1)) where NP = new{NP}(ϵ)
    SimplexCrossover{NP}(params::Parameters) where NP = new{NP}(params[:SPX_ϵ])
end

const SPX_DefaultOptions = ParamsDict(
    :SPX_ϵ => 1.1 # how much to inflate
)

#masscenter(pop, parentIndices) = mapslices(mean, pop[:, parentIndices], 2)

function apply!(xover::SimplexCrossover{NP},
                target::Individual, targetIndex::Int,
                pop, parentIndices) where NP
    @assert length(parentIndices) == NP

    offset = fill!(similar(target), 0)
    parentSum = copy(viewer(pop, parentIndices[1]))
    for i in 1:(NP-1)
        s = rand()^(1/i)
        parent = viewer(pop, parentIndices[i])
        parentNext = viewer(pop, parentIndices[i+1])
        @inbounds for j in 1:numdims(pop)
            offset[j] = s * (xover.ϵ*(parent[j] - parentNext[j]) + offset[j])
            parentSum[j] += parentNext[j]
        end
    end
    s = (1-xover.ϵ)/NP
    parentLast = viewer(pop, parentIndices[end])
    @inbounds for j in 1:numdims(pop)
        target[j] = s*parentSum[j] + xover.ϵ*parentLast[j] + offset[j]
    end
    return target
end
