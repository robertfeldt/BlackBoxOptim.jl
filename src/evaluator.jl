abstract Evaluator

type ProblemEvaluator <: Evaluator
  problem::OptimizationProblem
  fitness_scheme::FitnessScheme
  archive::Archive
  num_evals::Int
  last_fitness

  ProblemEvaluator(problem::OptimizationProblem; archive = false,
    fitness_scheme = ScalarFitness{true}()) = begin
    archive = archive || TopListArchive(numdims(problem))
    new(problem, fitness_scheme, archive, 0, nothing)
  end
end

worst_fitness(e::Evaluator) = worst_fitness(e.fitness_scheme)
num_evals(e::Evaluator) = e.num_evals
numdims(e::Evaluator) = numdims(e.problem)
search_space(e::Evaluator) = search_space(e.problem)
describe(e::Evaluator) = "Problem: $(name(e.problem)) (dimensions = $(numdims(e)))"
problem_summary(e::Evaluator) = "$(name(e.problem))_$(numdims(e))d"
problem(e::Evaluator) = e.problem

function evaluate(e::Evaluator, candidate)
  e.last_fitness = fitness(candidate, e.problem)
  e.num_evals += 1
  add_candidate!(e.archive, e.last_fitness, candidate, e.num_evals)
  e.last_fitness
end

# A way to get the fitness of the last evaluated candidate. Leads to nicer
# code if we can first test if better than or worse than existing candidates
# and only want the fitness itself if the criteria fulfilled.
last_fitness(e::Evaluator) = e.last_fitness

is_better{F}(f1::F, f2::F, e::Evaluator) = is_better(f1, f2, e.fitness_scheme)

is_better(candidate, fitness, e::Evaluator) = is_better(evaluate(e, candidate), fitness, e.fitness_scheme)

function best_of(candidate1, candidate2, e::Evaluator)
  f1 = evaluate(e, candidate1)
  f2 = evaluate(e, candidate2)
  if is_better(f1, f2, e.fitness_scheme)
    return candidate1, f1
  else
    return candidate2, f2
  end
end

# Candidate for the introduction into the population
type Candidate{F}
    params::Individual
    index::Int           # index of individual in the population, -1 if unassigned
    fitness::F           # fitness

    Candidate(params::Individual, index::Int = -1,
              fitness::F = NaN) = new(params, index, fitness)
end

function rank_by_fitness!{F}(e::Evaluator, candidates::Vector{Candidate{F}})
  # Note that we re-evaluate all candidates here. This might be wasteful and
  # we should cache if evaluations are costly.
  for i in eachindex(candidates)
    candidates[i].fitness = evaluate(e, candidates[i].params)
  end
  sort!(candidates; by = c -> c.fitness)
end

Base.copy{F}(c::Candidate{F}) = Candidate{F}(copy(c.params), c.index, c.fitness)

fitness_is_within_ftol(e::Evaluator, atol::Float64) = fitness_is_within_ftol(problem(e), best_fitness(e.archive), atol)
