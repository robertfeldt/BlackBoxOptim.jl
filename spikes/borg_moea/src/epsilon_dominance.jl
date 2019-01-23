function epsilon_dominates_clear(u, v, epsilon)
  veps = v .+ epsilon
  all(u .<= veps) && any(u .< veps)
end

function epsilon_dominates_fast(u, v, epsilon)
  pareto_dominates_fast(u, v .+ epsilon)
end

function pareto_dominates_clear(u, v)
  all(u .<= v) && any(u .< v)
end

function pareto_dominates_fast(u, v)
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

# Returns two booleans where the first indicates e-box dominance and the 2nd indicates
# epsilon progress
function epsilon_box_dominates_and_epsilon_progress(u, v, epsilon)
  uindex = floor(u ./ epsilon)
  vindex = floor(v ./ epsilon)
  inner_epsilon_box_dominates_and_epsilon_progress(
    u, uindex, (epsilon .* uindex), 
    v, vindex, (epsilon .* vindex))
end

epsilon_box_dominates(u, v, epsilon) = epsilon_box_dominates_and_epsilon_progress(u, v, epsilon)[1]
epsilon_box_progress(u, v, epsilon) = epsilon_box_dominates_and_epsilon_progress(u, v, epsilon)[2]

# Inner boolean expression when the intermediate results have already been calculated.
# This is to speed up processing when there is a whole archive of solutions to compare
# a new solution to.
function inner_epsilon_box_dominates_and_epsilon_progress(
  u, uindex, uindextimesepsilon, 
  v, vindex, vindextimesepsilon)

  if pareto_dominates_fast(uindex, vindex)
    return (true, true)
  elseif all(uindex .== vindex)
    smaller_distance = norm(u .- uindextimesepsilon) < norm(v .- vindextimesepsilon)
    return (smaller_distance, false)
  else
    return (false, false)
  end

end

# An EpsilonBoxArchive keeps a epsilon Pareto set of the best, non-dominated solutions
# found so far for a multi-objective search/optimization problem.
mutable struct EpsilonBoxArchive
  epsilon
  numsolutions::Int
  solutions::Matrix{Float64}       # solutions in this archive, one solution per column

  # For speedup we cache the boxes and leftpoints of all solutions. They have the same index as in
  # the solutions above.
  boxes::Matrix{Float64}      # the epsilon-boxes corresponding to each solution (same index as solutions)
  leftpoints::Matrix{Float64} # Position of "lower left" point of box (same index as solutions)

  function EpsilonBoxArchive(dimension::Integer, epsilon = 0.01, size = 100)
    new(
      epsilon,
      0,
      zeros(Float64, dimension, size),
      zeros(Float64, dimension, size),
      zeros(Float64, dimension, size)
      )
  end
end

add_if_dominates!(a::EpsilonBoxArchive)