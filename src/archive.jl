# An archive saves information about interesting points during an optimization
# run. It also saves (at least) the two last best fitness values since they
# can be used to estimate a confidence interval for how much improvement
# potential there is.
abstract Archive

# A top list archive saves a top list of the best performing (best fitness)
# candidates/individuals seen.
type TopListArchive <: Archive
  size::Int64           # Max size of top lists
  count::Int64          # Current number of values in the top lists

  fitnesses::Array{Float64,1}  # Top fitness values
  candidates::Array{Any, 1}    # Top candidates corresponding to top fitness values

  best_fitness_history::Array{Float64,1}

  numdims::Int64        # Number of dimensions in opt problem. Needed for confidence interval estimation.

  TopListArchive(numdims, size::Int64 = 10) = begin
    new(size, 0, Float64[], Any[], Float64[], int(numdims))
  end
end

best_candidate(a::Archive) = a.candidates[1]
best_fitness(a::Archive) = a.fitnesses[1]
last_top_fitness(a::Archive) = a.fitnesses[a.count]

should_enter_toplist(fitness::Float64, a::Archive) = fitness < last_top_fitness(a)

# Add a candidate with a fitness to the archive (if it is good enough).
function add_candidate!(a::TopListArchive, fitness::Float64, candidate::Array{Float64, 1})
  if a.count < a.size

    if a.count == 0 || fitness < best_fitness(a)
      push_to_fitness_history!(a, fitness)
    end

    push_then_sort_by_fitness!(fitness, candidate, a)
    a.count += 1

  elseif should_enter_toplist(fitness, a)

    if fitness < best_fitness(a)
      push_to_fitness_history!(a, fitness)
    end

    push_then_sort_by_fitness!(fitness, candidate, a)

    # Cut so only top list values are saved
    a.fitnesses = a.fitnesses[1:a.size]
    a.candidates = a.candidates[1:a.size]

  end
end

function push_to_fitness_history!(a::Archive, fitness)
  push!(a.best_fitness_history, fitness)
end

function push_then_sort_by_fitness!(fitness::Float64, candidate::Array{Float64, 1}, a::Archive)
  push!(a.fitnesses, fitness)
  push!(a.candidates, candidate)
  order = sortperm(a.fitnesses)
  a.fitnesses = a.fitnesses[order]
  a.candidates = a.candidates[order]
end

# Calculate the width of the confidence interval at a certain p-value.
# This is based on the paper:
#  Carvalho (2011), "Confidence intervals for the minimum of a
#    function using extreme value statistics"
#
# This means that the current estimate of the confidence interval for the minimum
# of the optimized function lies within the interval
#
#     ] l1 - w, l1 [
#
# with probability (1-p) as the number of sampled function points goes to infinity,
# where
#   w = width_of_confidence_interval(a, p)
#   l1 = best_fitness(a)
#
function width_of_confidence_interval(a::Archive, p = 0.05)
  if length(a.best_fitness_history) < 2
    return Inf
  else
    l1 = a.best_fitness_history[end]
    l2 = a.best_fitness_history[end-1]
    # We use abs below so it works also for maximization.
    abs(l2 - l1) / ( (1-p)^(-2/a.numdims) - 1)
  end
end

# We define the improvement potential as the percentage improvement 
# that can be expected from the current fitness value at a given p value.
# In theory, an optimization run should be terminated when this value is very 
# small, i.e. there is little improvement potential left in the run.
function fitness_improvement_potential(a::Archive, p = 0.01)
  width_of_confidence_interval(a, p) / abs(best_fitness(a))
end