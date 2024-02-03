module BlackBoxOptim

using Distributions, StatsBase, Random, LinearAlgebra, Printf, Distributed, Compat
using SpatialIndexing
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
        AbstractAsyncEvaluator, AbstractFitnessEvaluationJob,
        update_fitness!, async_update_fitness!, sync_update_fitness,

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
        ParamBounds, Individual, SearchSpace, FixedDimensionSearchSpace,
        RectSearchSpace, ContinuousRectSearchSpace, MixedPrecisionRectSearchSpace,
        numdims, dimmin, dimmax, dimdelta, dimrange, dimdigits,
        rand_individual, rand_individuals,

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

if !isdefined(Base, :get_extension)
    using Requires
end

module Utils
    using Random

    include("utilities/latin_hypercube_sampling.jl")
    include("utilities/assign_ranks.jl")
    include("utilities/halton_sequence.jl")
end

const SI = SpatialIndexing

include("search_space.jl")
include("parameters.jl")
include("fitness.jl")
include("ntuple_fitness.jl")
include("sliding_bitset.jl")

include("problem.jl")

include("frequency_adaptation.jl")

include("fit_individual.jl")
include("archive.jl")
include("archives/dominance_cone.jl")
include("archives/epsbox_archive.jl")

include("genetic_operators/genetic_operator.jl")

include("evaluator.jl")
include("parallel_evaluator.jl")
include("multithread_evaluator.jl")

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

# Quality Diversity optimization algorithms
#include("behavioral_space.jl")

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

# GUIs and front-ends (to really use it, one needs HTTP to enable BlackBoxOptimRealtimePlotServerExt)
include(joinpath("gui", "realtime_plot.jl"))

@static if !isdefined(Base, :get_extension)
    function __init__()
        @require Sockets = "6462fe0b-24de-5631-8697-dd941f90decc" begin
            @require HTTP = "cd3eb016-35fb-5094-929b-558a96fad6f3" begin
                try
                    include("../ext/BlackBoxOptimRealtimePlotServerExt.jl")
                catch err
                    println("Error during pre-compilation, when loading the gui extension, ", err)
                end
            end
        end
    end
end

end # module BlackBoxOptim
