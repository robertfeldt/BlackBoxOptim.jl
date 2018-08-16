module BlackBoxOptim

using Distributions, StatsBase, Random, LinearAlgebra, Printf, Distributed, Compat
using Printf: @printf, @sprintf
using Compat: String, view

export  Optimizer, AskTellOptimizer, SteppingOptimizer, PopulationOptimizer,
        bboptimize, bbsetup, compare_optimizers,

        DiffEvoOpt, de_rand_1_bin, de_rand_1_bin_radiuslimited,
        adaptive_de_rand_1_bin, adaptive_de_rand_1_bin_radiuslimited,

        SeparableNESOpt, separable_nes,
        XNESOpt, xnes, dxnes,

        # Parameters
        DictChain, Parameters, ParamsDictChain, ParamsDict,
        chain, flatten,

        # Fitness
        FitnessScheme,
        ScalarFitnessScheme, ComplexFitnessScheme,
        MinimizingFitnessScheme, MaximizingFitnessScheme,
        IndexedTupleFitness, TupleFitnessScheme, ParetoFitnessScheme,
        EpsDominanceFitnessScheme, EpsBoxDominanceFitnessScheme,
        fitness_type, fitness_eltype, numobjectives,
        is_minimizing, nafitness, isnafitness,
        hat_compare, is_better, is_worse, same_fitness,
        aggregate,

        # Evaluator
        #ProblemEvaluator,

        # Problems
        Problems,
        OptimizationProblem, FunctionBasedProblem,
        minimization_problem,
        name, fitness_scheme_type, fitness_scheme, search_space, numdims, opt_value,
        fitness_is_within_ftol, objfunc, fitness,

        # Problem factory/family
        FunctionBasedProblemFamily, MinimizationProblemFamily,
        fixed_dim_problem, instantiate,

        save_fitness_history_to_csv_file,

        # Archive
        TopListArchive, EpsBoxArchive,
        best_fitness, best_candidate,
        last_top_fitness, delta_fitness, capacity,
        width_of_confidence_interval, fitness_improvement_potential,

        # OptimizationResults
        minimum, f_minimum, iteration_converged, parameters, population, pareto_frontier, params,
        archived_fitness,

        # OptController
        numruns, lastrun, problem,

        # Search spaces
        ParamBounds, Individual, SearchSpace, FixedDimensionSearchSpace, ContinuousSearchSpace,
        RangePerDimSearchSpace, symmetric_search_space,
        numdims, mins, maxs, deltas, ranges, range_for_dim, diameters,
        rand_individual, rand_individuals, rand_individuals_lhs,

        # Population
        FitPopulation,
        popsize,

        # Genetic operators
        GeneticOperator, MutationOperator, CrossoverOperator, EmbeddingOperator,
        NoMutation, MutationClock, GibbsMutationOperator, UniformMutation,
        PolynomialMutation,
        FixedGeneticOperatorsMixture, FAGeneticOperatorsMixture,
        RandomBound,
        SimpleSelector, RadiusLimitedSelector,
        SimulatedBinaryCrossover, SimplexCrossover, UnimodalNormalDistributionCrossover,
        ParentCentricCrossover,

        numparents, numchildren, apply!, adjust!,

        # Utilities
        FrequencyAdapter, update!, frequencies,
        name

module Utils
    using Random

    include("utilities/latin_hypercube_sampling.jl")
    include("utilities/assign_ranks.jl")
    include("utilities/halton_sequence.jl")
end

include("search_space.jl")
include("parameters.jl")
include("fitness.jl")
include("ntuple_fitness.jl")

include("problem.jl")

include("frequency_adaptation.jl")

include("fit_individual.jl")
include("archive.jl")
include("archives/epsbox_archive.jl")

include("genetic_operators/genetic_operator.jl")

include("evaluator.jl")
include("parallel_evaluator.jl")

include("population.jl")
include("optimizer.jl")

# Different optimization algorithms/methods
include("random_search.jl")
include("differential_evolution.jl")
include("adaptive_differential_evolution.jl")
include("natural_evolution_strategies.jl")
include("dx_nes.jl")
include("resampling_memetic_search.jl")
include("simultaneous_perturbation_stochastic_approximation.jl")
include("generating_set_search.jl")
include("direct_search_with_probabilistic_descent.jl")

# multi-objective optimization algorithms
include("borg_moea.jl")

# End-user/top-level interface functions
include(joinpath("problems", "problem_family.jl"))
include("optimization_methods.jl")
include("default_parameters.jl")
include("optimization_result.jl")
include("opt_controller.jl")
include("bboptimize.jl")
include("compare_optimizers.jl")

# Problems for testing
include(joinpath("problems", "single_objective.jl"))
include(joinpath("problems", "multi_objective.jl"))

end # module BlackBoxOptim
