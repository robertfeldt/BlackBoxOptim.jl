"""
    Parent Centric Crossover (PCX).

    See
        Deb, K., Anand, A., and Joshi, D., "A Computationally Efficient Evolutionary Algorithm for Real-Parameter Optimization," Evolutionary Computation, vol. 10, no. 4, pp. 371-395, 2002.
"""
immutable ParentCentricCrossover{NP} <: CrossoverOperator{NP,1}
	ζ::Float64 # sd for the orthogonal directions defined by 2nd,3rd etc parents
    η::Float64 # sd for the direction of the 1st parent

    function ParentCentricCrossover(ζ::Number, η::Number)
        ζ > 0 || throw(ArgumentError("ζ must be positive"))
        η > 0 || throw(ArgumentError("η must be non-negative"))
        new(ζ, η)
    end
    ParentCentricCrossover(params::Parameters) = new(params[:PCX_ζ], params[:PCX_η])
end

const PCX_DefaultOptions = ParamsDict(
  :PCX_ζ => 0.1,
  :PCX_η => 0.5
)

function apply!{NP}(xover::ParentCentricCrossover{NP}, target::Individual, targetIndex::Int, pop, parentIndices)
    @assert length(parentIndices) == NP

    parents_centered = pop[:, parentIndices]
    center = mean(parents_centered, 2)
    broadcast!(-, parents_centered, parents_centered, center)
    # project other parents vectors orthogonal to
    # the subspace orthogonal to the selected parent
    tmp_mtx = A_mul_Bt(view(parents_centered, :, 1), view(parents_centered, :, 1))
    @inbounds for i in 1:size(tmp_mtx, 1)
        tmp_mtx[i,i] -= 1.0
    end
    other_parents_centered = tmp_mtx *
            view(parents_centered, :, 2:length(parentIndices))
    sd = mean(map(sqrt, sumabs2(other_parents_centered, 2)))
    if sd > 1E-8
        svdf = svdfact(other_parents_centered)
        svals_norm = norm(svdf.S)
        svals = svdf.S * sd * xover.ζ / svals_norm
        svals .*= rand(Normal(), length(svals))

        A_mul_B!(target, svdf.U, svals)
    else # degenerated
        fill!(target, zero(eltype(target)))
    end

    selParent = viewer(pop, parentIndices[1])
    selScale = rand(Normal()) * xover.η
    selParentOffset = view(parents_centered, :, 1)
    @inbounds @simd for i in eachindex(target)
        target[i] += selParent[i] + selScale * selParentOffset[i]
    end

    return target
end
