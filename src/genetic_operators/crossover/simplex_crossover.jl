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

function apply!(xover::SimplexCrossover,
                target::Individual, targetIndex::Int,
                pop, parentIndices)
    @assert length(parentIndices) == numparents(xover)

    offset = fill!(similar(target), 0)
    parentSum = copy(viewer(pop, parentIndices[1]))
    for i in 1:(numparents(xover)-1)
        s = rand()^(1/i)
        parent = viewer(pop, parentIndices[i])
        parentNext = viewer(pop, parentIndices[i+1])
        @inbounds offset .= s * (xover.ϵ * (parent .- parentNext) .+ offset)
        @inbounds parentSum .+= parentNext
    end
    s = (1-xover.ϵ)/numparents(xover)
    parentLast = viewer(pop, parentIndices[end])
    @inbounds target .= s .* parentSum + xover.ϵ .* parentLast .+ offset
    return target
end
