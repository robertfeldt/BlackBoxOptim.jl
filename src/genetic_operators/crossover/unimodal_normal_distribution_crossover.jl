"""
Unimodal Normal Distribution Crossover (UNDX).

See
    Kita, H., Ono, I., and Kobayashi, S., "Multi-parental Extension of the Unimodal Normal Distribution Crossover for Real-coded Genetic Algorithms," Proceedings of the 1999 Congress on Evolutionary Computation, pp. 1581-1588, 1999.
    Deb, K., Anand, A., and Joshi, D., "A Computationally Efficient Evolutionary Algorithm for Real-Parameter Optimization," Evolutionary Computation, vol. 10, no. 4, pp. 371-395, 2002.
"""
struct UnimodalNormalDistributionCrossover{NP} <: CrossoverOperator{NP,1}
	ζ::Float64 # sd for parent-defined orthogonal directions
    η::Float64 # sd for the directions orthogonal to the parent-defined ones, divided by √n

    function UnimodalNormalDistributionCrossover{NP}(ζ::Number, η::Number) where NP
        ζ > 0 || throw(ArgumentError("ζ must be positive"))
        η >= 0 || throw(ArgumentError("η must be non-negative"))
        new{NP}(ζ, η)
    end
    UnimodalNormalDistributionCrossover{NP}(params::Parameters) where NP =
        new{NP}(params[:UNDX_ζ], params[:UNDX_η])
end

const UNDX_DefaultOptions = ParamsDict(
    :UNDX_ζ => 0.5,
    :UNDX_η => 0.1
)

function apply!(xover::UnimodalNormalDistributionCrossover{NP},
                target::Individual, targetIndex::Int,
                pop, parentIndices) where NP
    @assert length(parentIndices) == NP

    parents = pop[:, parentIndices]
    center = mean(parents, dims=2)
    parents .-= center
    sd = mean(map(sqrt, sum(abs2, parents, dims=2)))
    svdf = svd(parents, full=false)
    svals = svdf.S / norm(svdf.S)
    @inbounds for (i,s) in enumerate(svals)
        svals[i] = sd*randn() * (
            abs(s) > 1E-8 ? # check if parent subspace; FIXME is the threshold correct?
            s * xover.ζ : xover.η )
    end

    mul!(target, svdf.U, svals)
    @inbounds @simd for i in eachindex(target)
        target[i] += center[i,1]
    end
    return target
end
