"""
Base class for archive-specific component
of the `OptimizationResults`.
"""
abstract type ArchiveOutput end

"""
Base class for method-specific component
of the `OptimizationResults`.
"""
abstract type MethodOutput end

struct DummyMethodOutput <: MethodOutput end

# no method-specific output by default
MethodOutput(method::Optimizer) = DummyMethodOutput()

"""
The results of running optimization method.

Returned by `run!(oc::OptRunController)`.
Should be compatible (on the API level) with the `Optim` package.
See `make_opt_results()`.
"""
mutable struct OptimizationResults
    method::String           # FIXME symbol instead or flexible?
    stop_reason::String      # FIXME turn into type hierarchy of immutable reasons with their attached info
    iterations::Int
    start_time::Float64           # time (seconds) optimization started
    elapsed_time::Float64         # time (seconds) optimization finished
    parameters::Parameters        # all user-specified parameters to bboptimize()
    f_calls::Int                  # total number of fitness function evaluations
    fit_scheme::FitnessScheme     # fitness scheme used by the archive
    archive_output::ArchiveOutput # archive-specific output
    method_output::MethodOutput   # method-specific output

    function OptimizationResults(ctrl, oc)
        new(
            string(oc.parameters[:Method]),
            stop_reason(ctrl),
            num_steps(ctrl),
            start_time(ctrl), elapsed_time(ctrl),
            oc.parameters,
            num_func_evals(ctrl),
            fitness_scheme(evaluator(ctrl).archive),
            ArchiveOutput(evaluator(ctrl).archive),
            MethodOutput(ctrl.optimizer))
    end
end

stop_reason(or::OptimizationResults) = or.stop_reason
iterations(or::OptimizationResults) = or.iterations
start_time(or::OptimizationResults) = or.start_time
elapsed_time(or::OptimizationResults) = or.elapsed_time
parameters(or::OptimizationResults) = or.parameters
f_calls(or::OptimizationResults) = or.f_calls

fitness_scheme(or::OptimizationResults) = or.fit_scheme
best_candidate(or::OptimizationResults) = or.archive_output.best_candidate
best_fitness(or::OptimizationResults) = or.archive_output.best_fitness
# FIXME doesn't work if there's no best candidate
numdims(or::OptimizationResults) = length(best_candidate(or))

function general_stop_reason(or::OptimizationResults)
    detailed_reason = stop_reason(or)

    if occursin(r"Fitness .* within tolerance .* of optimum", detailed_reason)
        return "Within fitness tolerance of optimum"
    end

    if occursin(r"Delta fitness .* below tolerance .*", detailed_reason)
        return "Delta fitness below tolerance"
    end

    return detailed_reason
end

# Alternative nomenclature that mimics Optim.jl more closely.
# FIXME should be it be enabled only for MinimizingFitnessScheme?
Base.minimum(or::OptimizationResults) = best_candidate(or)
f_minimum(or::OptimizationResults) = best_fitness(or)
# FIXME lookup stop_reason
iteration_converged(or::OptimizationResults) = iterations(or) >= parameters(or)[:MaxSteps]
isinterrupted(or::OptimizationResults) = stop_reason(or) == "InterruptException()"

"""
`TopListArchive`-specific components of the optimization results.
"""
struct TopListArchiveOutput{F,C} <: ArchiveOutput
    best_fitness::F
    best_candidate::C

    TopListArchiveOutput(archive::TopListArchive{F}) where F =
        new{F,Individual}(best_fitness(archive), best_candidate(archive))
end

ArchiveOutput(archive::TopListArchive) = TopListArchiveOutput(archive)

"""
Wrapper for `FrontierIndividual` that allows easy access to the problem fitness.
"""
struct FrontierIndividualWrapper{F,FA} <: FitIndividual{F}
    inner::FrontierIndividual{FA}
    fitness::F

    FrontierIndividualWrapper{F}(
            indi::FrontierIndividual{FA},
            fit_scheme::FitnessScheme{FA}) where {F, FA} =
        new{F, FA}(indi, convert(F, fitness(indi), fit_scheme))
end

params(indi::FrontierIndividualWrapper) = params(indi.inner)
archived_fitness(indi::FrontierIndividualWrapper) = fitness(indi.inner)

"""
`EpsBoxArchive`-specific components of the optimization results.
"""
struct EpsBoxArchiveOutput{N,F,FS<:EpsBoxDominanceFitnessScheme} <: ArchiveOutput
    best_fitness::NTuple{N,F}
    best_candidate::Individual
    frontier::Vector{FrontierIndividualWrapper{NTuple{N,F},IndexedTupleFitness{N,F}}} # inferred Pareto frontier
    fit_scheme::FS

    function EpsBoxArchiveOutput(archive::EpsBoxArchive{N,F}) where {N,F}
        fit_scheme = fitness_scheme(archive)
        new{N,F,typeof(fit_scheme)}(convert(NTuple{N,F}, best_fitness(archive), fit_scheme),
                                    best_candidate(archive),
                                    FrontierIndividualWrapper{NTuple{N,F}, IndexedTupleFitness{N,F}}[FrontierIndividualWrapper{NTuple{N,F}}(archive.frontier[i], fit_scheme)
                                                                                                     for i in findall(archive.frontier_isoccupied)],
                                    fit_scheme)
    end
end

ArchiveOutput(archive::EpsBoxArchive) = EpsBoxArchiveOutput(archive)

pareto_frontier(or::OptimizationResults) = or.archive_output.frontier

"""
`PopulationOptimizer`-specific components of the `OptimizationResults`.
Stores the final population.
"""
struct PopulationOptimizerOutput{P} <: MethodOutput
    population::P

    # FIXME PO is only required so that julia doesn't think it's an old inner ctor syntax
    PopulationOptimizerOutput(method::PO) where {PO<:PopulationOptimizer} =
        new{typeof(population(method))}(population(method))
end

MethodOutput(optimizer::PopulationOptimizer) = PopulationOptimizerOutput(optimizer)

population(or::OptimizationResults) = or.method_output.population
