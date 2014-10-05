function epsilon_dominates(u, v, epsilon)
  veps = v .+ epsilon
  all(u .<= veps) && any(u .< veps)
end

function pareto_dominates(u, v)
  all(u .<= v) && any(u .< v)
end