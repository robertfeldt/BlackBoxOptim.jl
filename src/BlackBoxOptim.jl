module BlackBoxOptim

using Distributions, StatsBase, Compat

export  Optimizer, AskTellOptimizer, SteppingOptimizer, PopulationOptimizer,
        bboptimize, compare_optimizers,

        DiffEvoOpt, de_rand_1_bin, de_rand_1_bin_radiuslimited,

        AdaptConstantsDiffEvoOpt, adaptive_de_rand_1_bin, adaptive_de_rand_1_bin_radiuslimited,

        SeparableNESOpt, separable_nes,
        XNESOpt, xnes,

        # Parameters
        Parameters, mergeparam,

        # Fitness
        FitnessScheme,
        ScalarFitnessScheme, ComplexFitnessScheme, VectorFitnessScheme,
        MinimizingFitnessScheme, MaximizingFitnessScheme,
        fitness_type, numobjectives,
        is_minimizing, nafitness, isnafitness,
        hat_compare, is_better, is_worse, same_fitness,
        vector_fitness_scheme_min, vector_fitness_scheme_max,

        # Evaluator
        #ProblemEvaluator,

        # Problems
        Problems,
        OptimizationProblem, FunctionBasedProblem,
        name, fitness_scheme, search_space, numdims, opt_value,
        fitness_is_within_ftol, objfunc, fitness,

        # Problem factory/family
        FunctionBasedProblemFamily, MinimizationProblemFamily,
        fixed_dim_problem,

        save_fitness_history_to_csv_file,

        # Archive
        TopListArchive, best_fitness, best_candidate,
        last_top_fitness, delta_fitness, capacity,
        width_of_confidence_interval, fitness_improvement_potential,

        # Search spaces
        ParamBounds, Individual, SearchSpace, FixedDimensionSearchSpace, ContinuousSearchSpace,
        RangePerDimSearchSpace, symmetric_search_space,
        numdims, mins, maxs, deltas, ranges, range_for_dim, diameters,
        rand_individual, rand_individuals, isinspace, rand_individuals_lhs,

        # Population
        FitPopulation,
        popsize,

        # Genetic operators
        GeneticOperator, MutationOperator, CrossoverOperator, EmbeddingOperator,
        NoMutation, MutationClock, GibbsMutationOperator, SimpleGibbsMutation, MutationMixture,
        RandomBound,
        SimpleSelector, RadiusLimitedSelector,

        name

# base abstract class for black-box optimization algorithms
abstract Optimizer

# SteppingOptimizer's do not have an ask and tell interface since they would be
# complex to implement if forced into that form.
abstract SteppingOptimizer <: Optimizer
evaluator(so::SteppingOptimizer) = so.evaluator

# optimizer using ask()/..eval fitness../tell!() sequence at each step
abstract AskTellOptimizer <: Optimizer

# population-based optimizers
abstract PopulationOptimizer <: AskTellOptimizer
population(popopt::PopulationOptimizer) = popopt.population
popsize(popopt::PopulationOptimizer) = popsize(population(popopt))

module Utils
  include("utilities/latin_hypercube_sampling.jl")
  include("utilities/assign_ranks.jl")
end

include("search_space.jl")
include("parameters.jl")
include("fitness.jl")

# Genetic Operators
include("genetic_operators/genetic_operator.jl")

include("frequency_adaptation.jl")
include("archive.jl")

include(joinpath("problems", "all_problems.jl"))
include(joinpath("problems", "problem_family.jl"))

include("evaluator.jl")

function setup(o::Optimizer, evaluator::Evaluator)
  # Do nothing, override if you need to setup prior to the optimization loop
end

function finalize(o::Optimizer, evaluator::Evaluator)
  # Do nothing, override if you need to finalize something after the optimization loop
end

# The standard name function converts the type of the optimizer to a string
# and strips off trailing "Opt".
function name(o::Optimizer)
  s = string(typeof(o))
  if s[end-2:end] == "Opt"
    return s[1:end-3]
  else
    return s
  end
end

# Population
include("population.jl")

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
include("random_search.jl")
include("differential_evolution.jl")
include("adaptive_differential_evolution.jl")
include("natural_evolution_strategies.jl")
include("resampling_memetic_search.jl")
include("simultaneous_perturbation_stochastic_approximation.jl")
include("generating_set_search.jl")
include("direct_search_with_probabilistic_descent.jl")

# End-user/interface functions
include("bboptimize.jl")

# Fitness
# include("fitness/fitness_types.jl") FIXME merge it with fitness.jl
include("fitness/pareto_dominance.jl")
include("fitness/epsilon_pareto_dominance.jl")
include("fitness/epsilon_box_dominance.jl")

# Problems for testing
include(joinpath("problems", "single_objective.jl"))

end # module BlackBoxOptim
