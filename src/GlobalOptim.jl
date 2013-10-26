module GlobalOptim

export  OptimizationProblem,
        Optimizer, PopulationOptimizer, 
        DEOpt,
        Problems,
        search_space, rand_population, optimize

abstract Optimizer
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
  tic()
  for(step in 1:numSteps)
    print(".")
    candidates = ask(opt)

    ranked_candidates = rank_by_fitness(candidates, problem)

    tell!(opt, ranked_candidates)
  end

  #show(opt)

  best, index, fitness = find_best_individual(problem, opt)
  print("\nBest candidate found: "); show(best)
  print("\nFitness: "); show(fitness)
  print("\n")
  t = toc()
  println("Steps per second = $(numSteps/t)")

  return best, fitness
end

end # module GlobalOptim