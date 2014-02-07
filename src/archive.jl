# An archive saves information about interesting points during an optimization
# run. It also saves (at least) the two last best fitness values since they
# can be used to estimate a confidence interval for how much improvement
# potential there is.
abstract Archive

# A top list archive saves a top list of the best performing (best fitness)
# candidates/individuals seen.
type TopListArchive <: Archive
  start_time::Float64   # Time when archive created, we use this to approximate the starting time for the opt...
  num_fitnesses::Int64  # Number of calls to add_candidate

  size::Int64           # Max size of top lists
  count::Int64          # Current number of values in the top lists

  fitnesses::Array{Float64,1}  # Top fitness values
  candidates::Array{Any, 1}    # Top candidates corresponding to top fitness values

  # We save a fitness history that we can later print to a csv file.
  # For each magnitude class (as defined by magnitude_class function below) we
  # we save the first entry of that class. The tuple saved for each magnitude
  # class is: (magnitude_class, time, num_fevals, fitness, width_of_confidence_interval)
  fitness_history::Array{((Float64, Float64), Float64, Int64, Float64, Float64),1}
  last_magnitude_class::(Float64, Float64)

  numdims::Int64        # Number of dimensions in opt problem. Needed for confidence interval estimation.

  TopListArchive(numdims, size::Int64 = 10) = begin
    new(time(), 0, size, 0, Float64[], Any[], 
      ((Float64, Float64), Float64, Int64, Float64, Float64)[], (0.0, 0.0),
      int(numdims))
  end
end

best_candidate(a::Archive) = a.candidates[1]
best_fitness(a::Archive) = a.fitnesses[1]
last_top_fitness(a::Archive) = a.fitnesses[a.count]

should_enter_toplist(fitness::Float64, a::Archive) = fitness < last_top_fitness(a)

# Delta fitness is the difference between the top two candidates found so far.
# 
function delta_fitness(a::Archive)
  if length(a.fitness_history) < 2
    Inf
  else
    abs(a.fitness_history[2][4] - a.fitness_history[1][4])
  end
end

# Add a candidate with a fitness to the archive (if it is good enough).
function add_candidate!(a::TopListArchive, fitness::Float64, 
  candidate, num_fevals::Int64 = -1)
  a.num_fitnesses += 1

  if a.count < a.size

    if a.count == 0 || fitness < best_fitness(a)
      push_to_fitness_history!(a, fitness, num_fevals)
    end

    push_then_sort_by_fitness!(fitness, candidate, a)
    a.count += 1

  elseif should_enter_toplist(fitness, a)

    if fitness < best_fitness(a)
      push_to_fitness_history!(a, fitness, num_fevals)
    end

    push_then_sort_by_fitness!(fitness, candidate, a)

    # Cut so only top list values are saved
    a.fitnesses = a.fitnesses[1:a.size]
    a.candidates = a.candidates[1:a.size]

  end
end

# The magnitude class of a value is a tuple indicating its sign and scale in a 
# tuple. This is used for filtering so that we only need to save one history 
# value per magnitude class instead of saving them all.
function magnitude_class(f::Float64)
  if f == 0.0
    (-1.0, 1e100)
  else
    (sign(f), floor(10.0*log10(abs(f)))/10.0)
  end
end

# Save fitness history so we can reconstruct the most important events later.
# We do this by only saving the first fitness event in its magnitude class, see
# above.
function push_to_fitness_history!(a::Archive, fitness::Float64, num_fevals::Int64 = -1)
  mc = magnitude_class(fitness)
  if mc != a.last_magnitude_class
    a.last_magnitude_class = mc
    nf = num_fevals == -1 ? a.num_fitnesses : num_fevals
    event = (mc, time(), nf, fitness, width_of_confidence_interval(a, 0.01))
    push!(a.fitness_history, event)
  end
end

function fitness_history_csv_header(a::Archive)
  "Date,Time,ElapsedTime,Magnitude,NumFuncEvals,ConfIntervalWidth0_01,Fitness"
end

function save_fitness_history_to_csv_file(a::Archive, filename = "fitness_history.csv";
  header_prefix = "", line_prefix = "", include_header = true)
  fh = open(filename, "a+")
  if include_header
    println(fh, join([header_prefix, fitness_history_csv_header(a)], ","))
  end
  for(event in a.fitness_history)
    mc, t, nf, f, w = event
    println(fh, join([line_prefix, strftime("%Y-%m-%d,%T", t), t-a.start_time,
      mc[2], nf, w, f], ","))
  end
  close(fh)
end

function push_then_sort_by_fitness!(fitness::Float64, candidate, a::Archive)
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
function width_of_confidence_interval(a::Archive, p = 0.01)
  if length(a.fitness_history) < 2
    return Inf
  else
    l1 = a.fitnesses[1]
    l2 = a.fitnesses[2]
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