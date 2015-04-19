function epsilon_dominates_clear(u, v, epsilon)
  veps = v .+ epsilon
  all(u .<= veps) && any(u .< veps)
end

function epsilon_dominates_fast(u, v, epsilon)
  pareto_dominates_fast(u, v .+ epsilon)
end
