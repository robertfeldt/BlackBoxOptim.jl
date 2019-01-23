pareto_dominates_clear(u::Vector{T}, v::Vector{T}) where {T <: Real} =
    all(u .<= v) && any(u .< v)

function pareto_dominates_fast(u::Vector{T}, v::Vector{T}) where {T <: Real}
    anylt = false
    @inbounds for i in eachindex(u)
        if u[i] > v[i]
            return false
        elseif u[i] < v[i]
            anylt = true
        end
    end
    return anylt
end

pareto_dominates(v1::Vector{T}, v2::Vector{T}) where {T <: Real} =
    pareto_dominates_fast(v1, v2)

function pareto_dominates_hat(u::Vector{T}, v::Vector{T}) where {T <: Real}
    anygt = anylt = false
    @inbounds for i in eachindex(u)
        if u[i] > v[i]
            anygt = true
        elseif u[i] < v[i]
            anylt = true
        end
    end
    if anygt
        return anylt ? 0 : 1
    else
        return anylt ? -1 : 0
    end
end

pareto_dominates(v1::Vector{T}, v2::Vector{T}) where {T <: Real} =
    pareto_dominates_fast(v1, v2)

# FIXME
pareto_dominates(f1::NewFitness, f2::NewFitness) =
    pareto_dominates(fitnessvalues(f1), fitnessvalues(f2))
