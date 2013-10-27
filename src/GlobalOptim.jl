module GlobalOptim

export  OptimizationProblem,

        Optimizer, PopulationOptimizer, 
        optimize,

        DiffEvoOpt, de_rand_1_bin, de_rand_1_bin_radiuslimited,

        AdaptConstantsDiffEvoOpt, adaptive_de_rand_1_bin, adaptive_de_rand_1_bin_radiuslimited,

        Problems,

        search_space, rand_population, latin_hypercube_sampling,
        hat_compare, isbetter, isworse, samefitness,
        popsize,
        FloatVectorFitness, float_vector_scheme_min, float_vector_scheme_max,
        FloatVectorPopulation

abstract Optimizer

include("fitness.jl")
include("population.jl")

abstract PopulationOptimizer <: Optimizer

# Our design is inspired by the object-oriented, ask-and-tell "optimizer API 
# format" as proposed in:
#
#  Collette, Y., N. Hansen, G. Pujol, D. Salazar Aponte and 
#  R. Le Riche (2010). On Object-Oriented Programming of Optimizers - 
#  Examples in Scilab. In P. Breitkopf and R. F. Coelho, eds.: 
#  Multidisciplinary Design Optimization in Computational Mechanics, Wiley, 
#  pp. 527-565.
#  https://www.lri.fr/~hansen/collette2010Chap14.pdf
#
# but since Julia is not OO this is more reflected in certain patterns of how
# to specify and call optimizers. The basic ask-and-tell pattern is:
#
#   while !optimizer.stop
#     x = ask(optimizer)
#     y = f(x)
#     optimizer = tell(optimizer, x, y)
#   end
#
# after which the best solutions can be found by:
#
#   yopt, xopt = best(optimizer)
#
# We have extended this paradigm with the use of an archive that saves 
# information on what we have learnt about the search space as well as the
# best solutions found. For most multi-objective optimization problems there
# is no single optimum. Instead there are many pareto optimal solutions.
# An archive collects information about the pareto optimal set or some 
# approximation of it. Different archival strategies can be implemented.

# Different optimization algorithms
include("differential_evolution.jl")
include("adaptive_differential_evolution.jl")

# Problems for testing
include(joinpath("problems", "all_problems.jl"))

function rand_population(populationSize, searchSpace::Array{(Float64, Float64)})
  dims = length(searchSpace)
  mins = [s[1] for s=searchSpace]
  maxs = [s[2] for s=searchSpace]
  deltas = maxs - mins
  # Basically min + delta * rand(), but broadcast over the columns...
  broadcast(+, mins', broadcast(*, deltas', rand(populationSize, dims)))
end

function find_best_individual(problem::Problems.OptimizationProblem, opt::PopulationOptimizer)
  pop = opt.population
  candidates = [(pop[i,:], i) for i in 1:size(pop,1)]
  rank_by_fitness(candidates, problem)[1]
end

function rank_by_fitness(candidates, problem)
  func = problem.funcs[1]
  # Note that we re-evaluate all candidates here. This might be wasteful and
  # we should cache if evaluations are costly.
  fitness = [(c[1], c[2], func(c[1])) for c=candidates]
  sort(fitness; by = (t) -> t[3])
end

function optimize(problem::Problems.OptimizationProblem, opt::Optimizer, numSteps = 1e4)
  num_better = 0
  println("----------------------------------------------------------------------")
  println("!!! Starting optimization !!!!")

  tic()
  for(step in 1:numSteps)
    if(mod(step, 2.5e4) == 0)
      println("Step $(step): Improvements/step = $(num_better/step)")
    end
    candidates = ask(opt)

    ranked_candidates = rank_by_fitness(candidates, problem)

    num_better += tell!(opt, ranked_candidates)
  end

  t = toc()
  println("Steps per second = $(numSteps/t)")

  if(mod(numSteps, 2.5e4) != 0)
    println("Step $(numSteps): Improvements/step = $(num_better/numSteps)")
  end

  println("\nMean value (in population) per position:")
  show(mean(opt.population,1))
  println("\n\nStd dev (in population) per position:")
  show(std(opt.population,1))

  best, index, fitness = find_best_individual(problem, opt)
  print("\n\nBest candidate found: "); show(best)
  print("\n\nFitness: "); show(fitness)
  println("\n----------------------------------------------------------------------")

  return best, fitness
end

end # module GlobalOptim