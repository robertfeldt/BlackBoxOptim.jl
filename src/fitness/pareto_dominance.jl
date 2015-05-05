pareto_dominates_clear{T <: Real}(u::Vector{T}, v::Vector{T}) = all(u .<= v) && any(u .< v)

function pareto_dominates_fast{T <: Real}(u::Vector{T}, v::Vector{T})
  one_smaller = false
  for i in 1:length(u)
    if u[i] > v[i]
      return false
    elseif u[i] < v[i]
      one_smaller = true
    end
  end
  return one_smaller
end

pareto_dominates{T <: Real}(v1::Vector{T}, v2::Vector{T}) = pareto_dominates_fast(v1, v2)

function pareto_dominates_hat{T <: Real}(u::Vector{T}, v::Vector{T})
  one_larger = one_smaller = false
  for i in 1:length(u)
    if u[i] > v[i]
      one_larger = true
    elseif u[i] < v[i]
      one_smaller = true
    end
  end
  if one_larger
    return (one_smaller ? 0 : 1)
  else
    return (one_smaller ? -1 : 0)
  end
end

pareto_dominates{T <: Real}(v1::Vector{T}, v2::Vector{T}) = pareto_dominates_fast(v1, v2)

pareto_dominates(f1::NewFitness, f2::NewFitness) = pareto_dominates(fitnessvalues(f1), fitnessvalues(f2))