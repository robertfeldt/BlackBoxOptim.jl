ValidMethods = {
  :random_search => BlackBoxOptim.random_search,
  :de_rand_1_bin => BlackBoxOptim.de_rand_1_bin,
  :de_rand_2_bin => BlackBoxOptim.de_rand_2_bin,
  :de_rand_1_bin_radiuslimited => BlackBoxOptim.de_rand_1_bin_radiuslimited,
  :de_rand_2_bin_radiuslimited => BlackBoxOptim.de_rand_2_bin_radiuslimited,
  :adaptive_de_rand_1_bin => BlackBoxOptim.adaptive_de_rand_1_bin,
  :adaptive_de_rand_1_bin_radiuslimited => BlackBoxOptim.adaptive_de_rand_1_bin_radiuslimited,
  :separable_nes => BlackBoxOptim.separable_nes,
  :xnes => BlackBoxOptim.xnes,
  :resampling_memetic_search => BlackBoxOptim.resampling_memetic_searcher,
  :resampling_inheritance_memetic_search => BlackBoxOptim.resampling_inheritance_memetic_searcher
}

function compare_optimizers(func::Function, searchRange; methods = keys(ValidMethods),
  iterations = 10000,
  max_time = 3.0,
  dimensions = :NotSpecified,
  show_trace::Bool = true,
  save_trace::Bool = false,
  population_size::Integer = 50,
  method_params = Dict())

  results = Any[]
  for(m in methods)
    tic()
    best, fitness = bboptimize(func, searchRange; method = m, iterations = iterations,
      max_time = max_time, dimensions = dimensions, 
      show_trace = show_trace, save_trace = save_trace, 
      population_size = population_size, method_params = method_params)
    push!( results,  (m, best, fitness, toq()) )
  end

  sorted = sort( results, by = (t) -> t[3] )

  if show_trace
    for(i in 1:length(sorted))
      println("$(sorted[i][1]), fitness = $(sorted[i][3]), time = $(sorted[i][4])")
    end
  end

  return sorted
end

function compare_optimizers(funcsAndRanges; methods = collect(keys(ValidMethods)),
  iterations = 10000,
  max_time = 2.0,
  dimensions = :NotSpecified,
  show_trace::Bool = true,
  save_trace::Bool = false,
  population_size::Integer = 50,
  method_params = Dict())
  
  # Lets create an array where we will save how the methods ranks per problem.
  ranks = zeros(length(methods), length(funcsAndRanges))
  fitnesses = zeros(Float64, length(methods), length(funcsAndRanges))

  for(i in 1:length(funcsAndRanges))
    res = compare_optimizers(funcsAndRanges[i][1], funcsAndRanges[i][2]; methods = methods, iterations = iterations,
      dimensions = dimensions, show_trace = show_trace, save_trace = save_trace, 
      population_size = population_size, method_params = method_params, max_time = max_time)
    for(j in 1:length(res))
      method, best, fitness, time = res[j]
      index = findfirst(methods, method)
      ranks[index, i] = j
      fitnesses[index, i] = fitness
    end
  end

  avg_ranks = mean(ranks, 2)
  avg_fitness = mean(fitnesses, 2)

  perm = sortperm(avg_ranks[:])
  println("By avg rank:")
  for(i in 1:length(methods))
    j = perm[i]
    print("$(i). $(methods[j]), average rank = $(avg_ranks[j]), avg fitness = $(avg_fitness[j]), ranks = $(ranks[j,:])")
  end

  perm = sortperm(avg_fitness[:])
  println("By avg fitness:")
  for(i in 1:length(methods))
    j = perm[i]
    print("$(i). $(methods[j]), average rank = $(avg_ranks[j]), avg fitness = $(avg_fitness[j]), ranks = $(ranks[j,:])")
  end

  return ranks, fitnesses
end

DefaultParameters = {
  :MaxTime => false,
  :MaxSteps => 10000,
  :Dimensions => :NotSpecified,
  :PopulationSize => 50,
  :FitnessTolerance => 1e-11
}

function bboptimize(func::Function, searchRange; method = :adaptive_de_rand_1_bin_radiuslimited,
  max_time = false,
  iterations = 10000,
  dimensions = :NotSpecified,
  show_trace::Bool = true,
  save_trace::Bool = false,
  population_size::Integer = 50,
  method_params = Dict())

  # Check that a valid search space has been stated and create the search_space
  # based on it, or bail out.
  if typeof(searchRange) == typeof( (0.0, 1.0) )
    if dimensions == :NotSpecified
      throw(ArgumentError("You MUST specify the number of dimensions in a solution when giving a search range $(searchRange)"))
    end
    search_space = symmetric_search_space(dimensions, searchRange)
  elseif typeof(searchRange) == typeof( [(0.0, 1.0)] )
    if dimensions != :NotSpecified
      throw(ArgumentError("You CANNOT specify the number of dimensions in a solution when first stating the search space ranges $(searchRange)"))
    end
    dimensions = length(searchRange)
    search_space = RangePerDimSearchSpace(searchRange)
  else
    throw(ArgumentError("Invalid search range specification. Either give only a range (ex: (0.0, 1.0)) AND the number of dimensions (ex: 2) or give an array of ranges (ex: [(0.0, 1.0), (10.0, 12.5)])."))
  end

  # Create a random solution from the search space and ensure that the given function returns a Number.
  ind = rand_individual(search_space)
  res = func(ind)
  if !(typeof(res) <: Number)
    throw(ArgumentError("The supplied function does NOT return a number when called with a potential solution (when called with $(ind) it returned $(res)) so we cannot optimize it!"))
  end

  # Check that max_time is larger than zero if it has been specified.
  if max_time != false
    if max_time <= 0.0
      throw(ArgumentError("The max_time must be a positive number"))
    else
      max_time = convert(Float64, max_time)
    end
  end

  # Check that a valid number of iterations has been specified. Print warning if higher than 1e8.
  if iterations < 1
    throw(ArgumentError("The number of iterations MUST be a positive number"))
  elseif iterations >= 1e8
    println("Number of allowed iterations is $(iterations); this can take a LONG time")
  end

  # Check that a valid population size has been given.
  if population_size < 2
    throw(ArgumentError("The population size MUST be at least 2"))
  end

  # Now create an optimization problem with the given information. We currently reuse the type
  # from our pre-defined problems so some of the data for the constructor is dummy.
  problem = BlackBoxOptim.FixedDimProblem("interactive", [func], search_space)

  # Check that a valid method has been specified and then set up the optimizer
  if (typeof(method) != Symbol) || !any([(method == vm) for vm in keys(ValidMethods)])
    throw(ArgumentError("The method specified, $(method), is NOT among the valid methods: $(ValidMethods)")) 
  end
  pop = BlackBoxOptim.rand_individuals_lhs(search_space, population_size)
  optimizer_func = ValidMethods[method]

  params = Parameters(method_params, DefaultParameters, {
    :Evaluator    => ProblemEvaluator(problem),
    :Population   => pop,
    :SearchSpace  => search_space
  })
  optimizer = optimizer_func(params)

  run_optimizer_on_problem(optimizer, problem; params = params,
   numSteps = iterations, shw = show_trace, save = save_trace, max_time = max_time)
end

function tr(msg, showTrace, saveTrace, obj = None)
  if showTrace
    print(msg)
    if obj != None
      show(obj)
    end
  end
  if saveTrace
    # No saving for now
  end
end

function find_best_individual(e::Evaluator, opt::PopulationOptimizer)
  (best_candidate(e.archive), 1, best_fitness(e.archive))
#  best_candidate(e.archive)
  #pop = population(opt)
  #candidates = [(pop[i,:], i) for i in 1:size(pop,1)]
  #rank_by_fitness(e, candidates)[1]
end

function find_best_individual(e::Evaluator, opt::Optimizer)
  (best_candidate(e.archive), 1, best_fitness(e.archive))
  #(opt.best, 1, opt.best_fitness)
end

function run_optimizer_on_problem(opt::Optimizer, problem::OptimizationProblem;
  numSteps = 1e4, 
  shw = true, 
  save = false, 
  max_time = false,
  fitness_tolerance = 1e-11,
  params = Dict())

  # No max time if unspecified. If max time specified it takes precedence over
  # numSteps.
  if max_time == false
    max_time = Inf
  else
    numSteps = Inf
  end

  # Now set up an evaluator for this problem. This will handle fitness
  # comparisons, keep track of the number of function evals as well as
  # keep an archive and top list.
  evaluator = get(params, :Evaluator, ProblemEvaluator(problem))

  num_better = 0
  num_better_since_last = 0
  tr("Starting optimization with optimizer $(name(opt))\n", shw, save)

  step = 1
  t = last_report_time = start_time = time()
  elapsed_time = 0.0

  while( true )

    if elapsed_time > max_time
      termination_reason = "Exceeded time budget"
      break
    end

    if step > numSteps
      termination_reason = "Exceeded max number of steps"
      break
    end

    if delta_fitness(evaluator.archive) < fitness_tolerance
      termination_reason = "Delta fitness below tolerance"
      break
    end

    # Report every 0.5 seconds
    if (t - last_report_time) > 0.5
      last_report_time = t
      num_better += num_better_since_last

      # Always print step number, num fevals and elapsed time
      tr(@sprintf("%.2f secs, %d evals , %d steps", 
        elapsed_time, num_evals(evaluator), step), shw, save) 

      # Only print if this optimizer reports on number of better. They return 0
      # if they do not.
      if num_better_since_last > 0
        tr(@sprintf(", improv/step: %.3f (last = %.4f)", 
          num_better/step, num_better_since_last/step), shw, save)
        num_better_since_last = 0
      end

      # Always print fitness if num_evals > 0
      if num_evals(evaluator) > 0
        tr(@sprintf(", %.9f", best_fitness(evaluator.archive)), shw, save)
      end

      tr("\n", shw, save)
    end

    if has_ask_tell_interface(opt)

      # They ask and tell interface is more general since you can mix and max
      # elements from several optimizers using it. However, in this top-level
      # execution function we do not make use of this flexibility...
      candidates = ask(opt)
      ranked_candidates = BlackBoxOptim.rank_by_fitness(evaluator, candidates)
      num_better_since_last += tell!(opt, ranked_candidates)

    else

      BlackBoxOptim.step(opt)
      num_better_since_last = 0

    end

    step += 1
    t = time()
    elapsed_time = t - start_time
  end

  step -= 1 # Since it is one too high after while loop above

  tr("\nOptimization stopped after $(step) steps and $(elapsed_time) seconds\n", shw, save)
  tr("Steps per second = $(step/elapsed_time)\n", shw, save)
  tr("Function evals per second = $(num_evals(evaluator)/elapsed_time)\n", shw, save)
  tr("Improvements/step = $((num_better+num_better_since_last)/numSteps)\n", shw, save)
  if typeof(opt) <: PopulationOptimizer
    tr("\nMean value (in population) per position:", shw, save, mean(population(opt),1))
    tr("\n\nStd dev (in population) per position:", shw, save, std(population(opt),1))
  end
  best, index, fitness = find_best_individual(evaluator, opt)
  tr("\n\nBest candidate found: ", shw, save, best)
  tr("\n\nFitness: ", shw, save, fitness)
  tr("\n\n", shw, false)

  if save
    optname = replace(name(optimizer), r"\s+", "_")
    archive = evaluator.archive
    # Save history to csv file here...
  end

  return best, fitness
end