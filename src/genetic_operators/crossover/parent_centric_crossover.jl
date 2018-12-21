"""
Parent Centric Crossover (PCX).

See
    Deb, K., Anand, A., and Joshi, D., "A Computationally Efficient Evolutionary Algorithm for Real-Parameter Optimization," Evolutionary Computation, vol. 10, no. 4, pp. 371-395, 2002.
"""
struct ParentCentricCrossover{NP} <: CrossoverOperator{NP,1}
    ζ::Float64 # sd for the orthogonal directions defined by 2nd,3rd etc parents
    η::Float64 # sd for the direction of the 1st parent

    function ParentCentricCrossover{NP}(ζ::Number, η::Number) where NP
        ζ > 0 || throw(ArgumentError("ζ must be positive"))
        η > 0 || throw(ArgumentError("η must be non-negative"))
        new{NP}(ζ, η)
    end
    ParentCentricCrossover{NP}(params::Parameters) where NP =
        new{NP}(params[:PCX_ζ], params[:PCX_η])
end

const PCX_DefaultOptions = ParamsDict(
    :PCX_ζ => 0.1,
    :PCX_η => 0.5
)

function apply!(xover::ParentCentricCrossover,
                target::Individual, targetIndex::Int,
                pop, parentIndices)
    @assert length(parentIndices) == numparents(xover)

    center = mean(viewer(pop, parentIndices), dims=2)
    # project other parents vectors to the subspace
    # orthogonal to the 1st (central) parent
    parent_1 = viewer(pop, parentIndices[1]) .- dropdims(center, dims=2)
    other_parents_centered = viewer(pop, parentIndices[2:end]) .- center
    if (parent_1_sqnorm = sum(abs2, parent_1)) > 0.0
        dotprods = transpose(parent_1) * other_parents_centered
        dotprods ./= parent_1_sqnorm
        other_parents_centered .-= dotprods .* parent_1
    end
    sd = mean(sqrt.(sum(abs2, other_parents_centered, dims=2)))
    if sd > 1E-8
        svdf = svd(other_parents_centered, full=false)
        svals_norm = norm(svdf.S)
        svals = svdf.S * sd * xover.ζ / svals_norm
        svals .*= randn(length(svals))
        mul!(target, svdf.U, svals)
    else # degenerated
        fill!(target, zero(eltype(target)))
    end

    selParent = viewer(pop, parentIndices[1])
    selScale = randn() * xover.η
    @inbounds target .+= selParent .+ selScale .* parent_1
    return target
end
