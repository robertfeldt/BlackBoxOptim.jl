"""
    Unimodal Normal Distribution Crossover (UNDX).

    See
        Kita, H., Ono, I., and Kobayashi, S., "Multi-parental Extension of the Unimodal Normal Distribution Crossover for Real-coded Genetic Algorithms," Proceedings of the 1999 Congress on Evolutionary Computation, pp. 1581-1588, 1999.
        Deb, K., Anand, A., and Joshi, D., "A Computationally Efficient Evolutionary Algorithm for Real-Parameter Optimization," Evolutionary Computation, vol. 10, no. 4, pp. 371-395, 2002.
"""
immutable UnimodalNormalDistributionCrossover{NP} <: CrossoverOperator{NP,1}
	ζ::Float64 # sd for parent-defined orthogonal directions
    η::Float64 # sd for the directions orthogonal to the parent-defined ones, divided by √n

    function UnimodalNormalDistributionCrossover(ζ::Number, η::Number)
        ζ > 0 || throw(ArgumentError("ζ must be positive"))
        η >= 0 || throw(ArgumentError("η must be non-negative"))
        new(ζ, η)
    end
    UnimodalNormalDistributionCrossover(params::Parameters) = new(params[:UNDX_ζ], params[:UNDX_η])
end

const UNDX_DefaultOptions = ParamsDict(
  :UNDX_ζ => 0.5,
  :UNDX_η => 0.1
)

function apply!{NP}(xover::UnimodalNormalDistributionCrossover{NP}, target::Individual, targetIndex::Int, pop, parentIndices)
    @assert length(parentIndices) == NP

    parents = pop[:, parentIndices]
    center = mean(parents, 2)
    broadcast!(-, parents, parents, center)
    sd = mean(map(sqrt, sumabs2(parents, 2)))
    svdf = svdfact(parents)
    svals = svdf.S / norm(svdf.S)
    parent_mask = Bool[abs(s) > 1E-8 for s in svals] # FIXME is the threshold correct?
    svals[parent_mask] .*= sd * xover.ζ
    svals[!parent_mask] = sd * xover.η
    svals .*= rand(Normal(), length(svals))

    A_mul_B!(target, svdf.U, svals)
    @inbounds @simd for i in eachindex(target)
        target[i] += center[i,1]
    end
    return target
end
