ValidMethods = {
  :random_search => BlackBoxOptim.random_search,
  :de_rand_1_bin => BlackBoxOptim.de_rand_1_bin,
  :de_rand_2_bin => BlackBoxOptim.de_rand_2_bin,
  :adaptive_de_rand_1_bin => BlackBoxOptim.adaptive_de_rand_1_bin,
  :de_rand_1_bin_radiuslimited => BlackBoxOptim.de_rand_1_bin_radiuslimited,
  :adaptive_de_rand_1_bin_radiuslimited => BlackBoxOptim.adaptive_de_rand_1_bin_radiuslimited,
  :separable_nes => BlackBoxOptim.separable_nes,
}

function compare_optimizers(func::Function, searchRange; methods = keys(ValidMethods),
  iterations = 10000,
  dimensions = :NotSpecified,
  show_trace::Bool = true,
  save_trace::Bool = false,
  population_size::Integer = 50,
  method_options = {})

  results = Any[]
  for(m in methods)
    tic()
    best, fitness = bboptimize(func, searchRange; method = m, iterations = iterations,
      dimensions = dimensions, show_trace = show_trace, save_trace = save_trace, 
      population_size = population_size, method_options = method_options)
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
  dimensions = :NotSpecified,
  show_trace::Bool = true,
  save_trace::Bool = false,
  population_size::Integer = 50,
  method_options = {})
  
  # Lets create an array where we will save how the methods ranks per problem.
  ranks = zeros(length(methods), length(funcsAndRanges))

  for(i in 1:length(funcsAndRanges))
    res = compare_optimizers(funcsAndRanges[i][1], funcsAndRanges[i][2]; methods = methods, iterations = iterations,
      dimensions = dimensions, show_trace = show_trace, save_trace = save_trace, 
      population_size = population_size, method_options = method_options)
    for(j in 1:length(res))
      method, best, fitness, time = res[j]
      index = findfirst(methods, method)
      ranks[index, i] = j
    end
  end

  avg_ranks = mean(ranks, 2)
  for(i in 1:length(methods))
    println("$(methods[i]), average rank = $(avg_ranks[i]), ranks = $(ranks[i,:])")
  end

  return ranks
end

function bboptimize(func::Function, searchRange; method = :adaptive_de_rand_1_bin_radiuslimited,
  iterations = 10000,
  dimensions = :NotSpecified,
  show_trace::Bool = true,
  save_trace::Bool = false,
  population_size::Integer = 50,
  method_options = {})

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

  # Check that a valid method has been specified and then set up the optimizer
  if (typeof(method) != Symbol) || !any([(method == vm) for vm in keys(ValidMethods)])
    throw(ArgumentError("The method specified, $(method), is NOT among the valid methods: $(ValidMethods)")) 
  end
  pop = BlackBoxOptim.rand_individuals_lhs(search_space, population_size)
  optimizer_func = ValidMethods[method]
  optimizer = optimizer_func(search_space; population = pop)

  # Now create an optimization problem with the given information. We currently reuse the type
  # from our pre-defined problems so some of the data for the constructor is dummy.
  problem = BlackBoxOptim.Problems.OptimizationProblem("interactive", [func], false, (0.0, 0.0), dimensions, search_space)

  run_optimizer_on_problem(optimizer, problem, iterations, show_trace, save_trace)
end

function tr(msg, showTrace, saveTrace, obj = None)
  if showTrace
    println(msg)
    if obj != None
      show(obj)
    end
  end
  if saveTrace
    # No saving for now
  end
end

function find_best_individual(problem::Problems.OptimizationProblem, opt::PopulationOptimizer)
  pop = population(opt)
  candidates = [(pop[i,:], i) for i in 1:size(pop,1)]
  rank_by_fitness(candidates, problem)[1]
end

function find_best_individual(problem::Problems.OptimizationProblem, opt::Optimizer)
  (opt.best, 1, opt.best_fitness)
end

function rank_by_fitness(candidates, problem)
  func = problem.funcs[1]
  # Note that we re-evaluate all candidates here. This might be wasteful and
  # we should cache if evaluations are costly.
  fitness = [(c[1], c[2], func(c[1])) for c=candidates]
  sort(fitness; by = (t) -> t[3])
end

function run_optimizer_on_problem(opt::Optimizer, problem::Problems.OptimizationProblem, 
  numSteps = 1e4, shw = true, save = false)

  num_better = 0
  num_better_since_last = 0
  tr("Starting optimization with optimizer $(opt.name)", shw, save)

  step = 1
  tic()
  while(step <= numSteps)
    if(mod(step, 2.5e4) == 0)
      num_better += num_better_since_last
      # Only print if this optimizer reports on number of better. They return 0
      # if they do not.
      if num_better_since_last > 0
        tr("Step $(step), Improvements/step: overall = $(num_better/step), last interval = $(num_better_since_last/step)", shw, save)
        num_better_since_last = 0
      end
    end
    candidates = ask(opt)

    ranked_candidates = rank_by_fitness(candidates, problem)

    num_better_since_last += tell!(opt, ranked_candidates)
    step += 1
  end
  t = toq()

  step -= 1 # Since it is one too high after while loop above

  tr("\nOptimization stopped after $(step) steps and $(t) seconds", shw, save)
  tr("Steps per second = $(numSteps/t)", shw, save)
  tr("Improvements/step = $((num_better+num_better_since_last)/numSteps)", shw, save)
  if typeof(opt) <: PopulationOptimizer
    tr("\nMean value (in population) per position:", shw, save, mean(population(opt),1))
    tr("\n\nStd dev (in population) per position:", shw, save, std(population(opt),1))
  end
  best, index, fitness = find_best_individual(problem, opt)
  tr("\n\nBest candidate found: ", shw, save, best)
  tr("\n\nFitness: ", shw, save, fitness)
  tr("\n", shw, false)

  return best, fitness
end