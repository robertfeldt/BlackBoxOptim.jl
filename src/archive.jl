# An archive saves information about interesting points during an optimization
# run. It also saves (at least) the two last best fitness values since they
# can be used to estimate a confidence interval for how much improvement
# potential there is.
abstract Archive

# A top list archive saves a top list of the best performing (best fitness)
# candidates/individuals seen.
type TopListArchive <: Archive
  start_time::Float64   # Time when archive created, we use this to approximate the starting time for the opt...
  num_fitnesses::Int  # Number of calls to add_candidate

  size::Int           # Max size of top lists
  count::Int          # Current number of values in the top lists

  fitnesses::Array{Float64,1}  # Top fitness values
  candidates::Array{Any, 1}    # Top candidates corresponding to top fitness values

  # We save a fitness history that we can later print to a csv file.
  # For each magnitude class (as defined by magnitude_class function below) we
  # we save the first entry of that class. The tuple saved for each magnitude
  # class is: (magnitude_class, time, num_fevals, fitness, width_of_confidence_interval)
  fitness_history::Array{(Float64, Int, Float64, Float64),1}

  numdims::Int        # Number of dimensions in opt problem. Needed for confidence interval estimation.

  TopListArchive(numdims, size::Int = 10) = begin
    new(time(), 0, size, 0, Float64[], Any[],
      (Float64, Int, Float64, Float64)[], int(numdims))
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
    abs(a.fitness_history[2][3] - a.fitness_history[1][3])
  end
end

# Add a candidate with a fitness to the archive (if it is good enough).
function add_candidate!(a::TopListArchive, fitness::Number,
  candidate, num_fevals::Int = -1)
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
function magnitude_class(f::Number)
  f = float(f)
  if f == 0.0
    (-1.0, 1e100)
  else
    (sign(f), floor(10.0*log10(abs(f)))/10.0)
  end
end

# Save fitness history so we can reconstruct the most important events later.
function push_to_fitness_history!(a::Archive, fitness::Number, num_fevals::Int = -1)
  push!(a.fitness_history, make_fitness_event(a, fitness, num_fevals))
end

function fitness_improvement_ratio(a::Archive, newFitness)
  try
    bestfitness = best_fitness(a)
    return abs( (bestfitness - newFitness) / bestfitness )
  catch
    return 0.0
  end
end

# Calculate the distance from a fitness value to the optimum/best known fitness value.
function distance_to_optimum(fitness, bestfitness)
  abs(fitness - bestfitness)
end

function make_fitness_event(a::Archive, fitness::Number, num_fevals::Int = -1)
  nf = num_fevals == -1 ? a.num_fitnesses : num_fevals
  (time(), nf, fitness, fitness_improvement_ratio(a, fitness))
end

function fitness_history_csv_header(a::Archive)
  "Date,Time,ElapsedTime,Magnitude,NumFuncEvals,FitnessImprovementRatio,Fitness"
end

function save_fitness_history_to_csv_file(a::Archive, filename = "fitness_history.csv";
  header_prefix = "", line_prefix = "",
  include_header = true, bestfitness = nothing)

  fh = open(filename, "a+")

  if include_header

    header = [header_prefix, fitness_history_csv_header(a)]

    if bestfitness != nothing
      push!(header, "DistanceToBestFitness")
    end

    println(fh, join(header, ","))

  end

  for(event in a.fitness_history)

    t, nf, f, fir = event

    mc = magnitude_class(f)

    line = [line_prefix, strftime("%Y-%m-%d,%T", t), t-a.start_time,
      mc[1]*mc[2], nf, fir, f]

    if bestfitness != nothing
      push!(line, distance_to_optimum(f, bestfitness))
    end

    println(fh, join(line, ","))

  end

  close(fh)

end

function push_then_sort_by_fitness!(fitness::Number, candidate, a::Archive)
  push!(a.fitnesses, fitness)
  push!(a.candidates, candidate)
  order = sortperm(a.fitnesses)
  a.fitnesses = a.fitnesses[order]
  a.candidates = a.candidates[order]
end

# Merge multiple fitness histories and calculate the min, max, avg and median
# values for fitness and fir at all fitness eval change points.
function merge_fitness_histories(histories)
  numhist = length(histories)
  counts = ones(numhist)
  fitnesses = zeros(numhist)
  firs = zeros(numhist)
  current_feval = 1

  # Find max history length
  histlengths = [length(h) for h in histories]
  maxlen = maximum(histlengths)

  while maximum(counts) < maxlen

    # Find min feval for current events, this is the next current feval
    for i in 1:numhist
      t, nf, f, fir = histories[i][counts[i]]
    end

  end

end

# Merge the fitness histories and save the average values of the fitness,
# and distance to best fitness for each change in any of the histories.
#function merge_fitness_histories_to_csv(archives::Archive[], filename = "fitness_history.csv";
#  header_prefix = "", line_prefix = "",
#  include_header = true, bestfitness = nothing)
#
#end

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
