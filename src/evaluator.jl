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

function update_fitness!(e::ProblemEvaluator{FP,FA}, candidate::Candidate{FA}) where {FP,FA}
    # evaluate fitness if not known yet
    if isnafitness(candidate.fitness, fitness_scheme(e.archive))
        candidate.fitness = fitness(candidate.params, e, candidate.tag)
    end
    return candidate
end

function update_fitness!(e::ProblemEvaluator{FP,FA}, candidates::AbstractVector{Candidate{FA}}) where {FP,FA}
    @inbounds for candidate in candidates
        update_fitness!(e, candidate)
    end
    return candidates
end

function rank_by_fitness!(e::Evaluator, candidates::AbstractVector{<:Candidate})
    fs = fitness_scheme(e)
    sort!(update_fitness!(e, candidates);
          # FIXME use lt=fitness_scheme(a) when v0.5 #14919 would be fixed
          by=fitness, lt=(x, y) -> is_better(x, y, fs))
end

# called by check_stop_condition(OptRunController)
check_stop_condition(e::Evaluator, ctrl) = check_stop_condition(e.archive, e.problem, ctrl)
