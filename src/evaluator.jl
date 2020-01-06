"""
The abstract base for types that manage the objective function evaluation.
`P` is the optimization problem it is used for.
"""
abstract type Evaluator{P <: OptimizationProblem} end

fitness_scheme(e::Evaluator) = fitness_scheme(problem(e))
fitness_type(e::Evaluator) = fitness_type(fitness_scheme(e))
worst_fitness(e::Evaluator) = worst_fitness(fitness_scheme(e))
numdims(e::Evaluator) = numdims(problem(e))
search_space(e::Evaluator) = search_space(problem(e))
describe(e::Evaluator) = "Problem: $(name(e.problem)) (dimensions = $(numdims(e)))"
problem_summary(e::Evaluator) = "$(name(e.problem))_$(numdims(e))d"

shutdown!(e::Evaluator) = e # do nothing

"""
    update_fitness!([f], eval::Evaluator, candidate::Candidate; force::Bool=false) -> Candidate

Calculate fitness of `candidate` and optionally apply `f`.
`force` specifies whether to re-evaluate fitness, if the candidate already has non-NA one.
"""
function update_fitness!(f::Any, e::Evaluator, candidate::Candidate; force::Bool=false)
    # evaluate fitness if not known yet
    if force || isnafitness(candidate.fitness, fitness_scheme(e.archive))
        candidate.fitness = fitness(candidate.params, e, candidate.tag)
        (f !== nothing) && f(candidate)
    end
    return candidate
end

update_fitness!(e::Evaluator, candidate::Candidate; force::Bool=false) =
    update_fitness!(nothing, e, candidate, force=force)

"""
    update_fitness!([f], eval::Evaluator, candidates; force::Bool=false)

Calculate fitness of `candidates` and optionally apply `f` to each processed one.
`force` specifies if already existing non-NA fitnesses should be re-evaluated.
"""
function update_fitness!(f::Any, e::Evaluator, candidates::Any; force::Bool=false)
    fs = fitness_scheme(e.archive)
    for candi in candidates
        if force || isnafitness(fitness(candi), fs)
            update_fitness!(f, e, candi)
        end
    end
    return candidates
end

update_fitness!(e::Evaluator, candidates::Any; force::Bool=false) =
    update_fitness!(nothing, e, candidates, force=force)

"""
The abstract base type for asynchronous evaluators.

Asynchronous evaluator allows submitting fitness calculation requests via
[`async_update_fitness!`](@ref) method that returns
[`AbstractFitnessEvaluationJob`](@ref) object without waiting for
the calculation completion. Job object could then be passed to
[`sync_update_fitness`](@ref) that waits until fitnesses for all specified jobs
is calculated.
"""
abstract type AbstractAsyncEvaluator{P <: OptimizationProblem} <: Evaluator{P} end

"""
The abstract base type for the fitness evaluation jobs returned by
[`async_update_fitness!`](@ref) and accepted by [`sync_update_fitness`](@ref).

One can also check whether the fitness completed with [`Base.isready`](@ref) method.
"""
abstract type AbstractFitnessEvaluationJob{F} end

"""
    async_update_fitness!(eval::AbstractAsyncEvaluator, candidates::Any;
                          force::Bool=false) ->
                          Union{AbstractFitnessEvaluationJob, Nothing}

Asynchronously calculate the fitnesses.
*candidates* should support the `iterate()` interface and return
[`Candidate`](@ref) elements.

Returns [`AbstractFitnessEvaluationJob`](@ref) or `nothing` if *candidates*
does not contain any candidates.

# Arguments
 - `force::Bool`: force fitness calculation even if the candidate already has
   non-NA fitness

See also: [`sync_update_fitness`](@ref)
"""
async_update_fitness!

"""
    sync_fitness_update([f], job::AbstractFitnessEvaluationJob,
                        eval::AbstractAsyncEvaluator)

Waits until fitnesses for `job` are calculated and optionally applies `f`
function to each processed candidate.

See also: [`async_update_fitness!`](@ref), [`update_fitness!`](@ref).
"""
sync_update_fitness

"""
    is_stopping(eval::AbstractAsyncEvaluator) -> Bool

Check if the evaluator is in the shutdown sequence.

When the evaluator is being shut down, it waits for the workers to finish
their jobs and stop, no new jobs could be submitted.
"""
is_stopping

"""
Default implementation of the `Evaluator`.

`FP` is the original problem's fitness type
`FA` is the fitness type actually stored by the archive.
"""
mutable struct ProblemEvaluator{FP, FA, A<:Archive, P<:OptimizationProblem} <: Evaluator{P}
    # FIXME F is the fitness type of the problem, but with current Julia it's
    # not possible to get and use it at declaration type
    problem::P
    archive::A
    num_evals::Int
    last_fitness::FP

    ProblemEvaluator(problem::P, archive::A) where {P<:OptimizationProblem, A<:Archive} =
        new{fitness_type(fitness_scheme(problem)),fitness_type(archive),A,P}(
                problem, archive,
                0, nafitness(fitness_scheme(problem)))
end

ProblemEvaluator(problem::OptimizationProblem; archiveCapacity::Int = 10) =
    ProblemEvaluator(problem, TopListArchive(fitness_scheme(problem), numdims(problem), archiveCapacity))

problem(e::Evaluator) = e.problem
num_evals(e::ProblemEvaluator) = e.num_evals

"""
    fitness(params::Individual, e::ProblemEvaluator, tag::Int=0)

Evaluate the fitness and implicitly update the archive with the provided
parameters and calculated fitness.

Returns the fitness in the archived format.
"""
function fitness(params::Individual, e::ProblemEvaluator, tag::Int=0)
    e.last_fitness = fit = fitness(params, e.problem)
    e.num_evals += 1
    fita = archived_fitness(fit, e.archive)
    candi = add_candidate!(e.archive, fita, params, tag, e.num_evals)
    return fita
end

"""
    last_fitness(e::Evaluator)

Get the fitness of the last evaluated candidate.

Leads to nicer code if we can first test if it is better or worse than existing candidates
and only want the fitness itself if the criteria fulfilled.
"""
last_fitness(e::Evaluator) = e.last_fitness

is_better(f1::F, f2::F, e::Evaluator) where {F} = is_better(f1, f2, fitness_scheme(e))

is_better(candidate, f, e::Evaluator) = is_better(fitness(candidate, e), f, fitness_scheme(e))

function best_of(candidate1::Individual, candidate2::Individual, e::Evaluator)
    f1 = fitness(candidate1, e)
    f2 = fitness(candidate2, e)
    if is_better(f1, f2, e)
        return candidate1, f1
    else
        return candidate2, f2
    end
end

function rank_by_fitness!(e::Evaluator, candidates::AbstractVector{<:Candidate})
    fs = fitness_scheme(e)
    sort!(update_fitness!(e, candidates);
          # FIXME use lt=fitness_scheme(a) when v0.5 #14919 would be fixed
          by=fitness, lt=(x, y) -> is_better(x, y, fs))
end

# called by check_stop_condition(OptRunController)
check_stop_condition(e::Evaluator, ctrl) = check_stop_condition(e.archive, e.problem, ctrl)
