module GlobalOptim

export  SingleObjectiveProblems, OptimizationProblem,
        Optimizer, PopulationOptimizer, 
        DEOpt

abstract Optimizer
abstract PopulationOptimizer <: Optimizer

# We base our design on the object-oriented, ask-and-tell "API format" for 
# writing optimizers as proposed in:
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

# Duplicate a tuple to indicate a whole search space, i.e. symmetrically.
function symmetric_search_space(dims, dimRange::(Float64, Float64))
  [dimRange for i=1:dims]
end

function rand_population(populationSize, searchSpace::Array{(Float64, Float64)})
  dims = length(searchSpace)
  mins = [s[1] for s=searchSpace]
  maxs = [s[2] for s=searchSpace]
  deltas = maxs - mins
  # Basically min + delta * rand(), but broadcast over the columns...
  broadcast(+, mins', broadcast(*, deltas', rand(populationSize, dims)))
end

function optimize(problem::Problems.OptimizationProblem, opt::Optimizer, numSteps = 1e4)
  for(step in 1:numSteps)
    candidates = ask(de)
    fitness = [(problem.f(c[1]), c[1], c[2]) for c=candidates]
    sorted = sort(fitness; by = (t) -> t[1])
    ranked_candidates = map((t) -> (t[2], t[3]), sorted)
    tell(de, ranked_candidates)
  end
end

end # module GlobalOptim