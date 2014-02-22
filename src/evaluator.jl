abstract Evaluator

type ProblemEvaluator <: Evaluator
  problem::OptimizationProblem
  fitness_scheme::FitnessScheme
  archive::Archive
  num_evals::Int64
  last_fitness

  ProblemEvaluator(problem::OptimizationProblem; archive = false,
    fitness_scheme = FloatFitness()) = begin
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
  e.last_fitness = eval1(candidate, e.problem)
  e.num_evals += 1
  add_candidate!(e.archive, e.last_fitness, candidate, e.num_evals)
  e.last_fitness
end

# A way to get the fitness of the last evaluated candidate. Leads to nicer
# code if we can first test if better than or worse than existing candidates
# and only want the fitness itself if the criteria fulfilled.
last_fitness(e::Evaluator) = e.last_fitness

function is_better(e::Evaluator, f1::Fitness, f2::Fitness)
  is_better(f1, f2, e.fitness_scheme)
end

function is_better(e::Evaluator, candidate, fitness::Fitness)
  is_better(evaluate(e, candidate), fitness, e.fitness_scheme)
end

function best_of(e::Evaluator, candidate1, candidate2)
  f1 = evaluate(e, candidate1)
  f2 = evaluate(e, candidate2)
  if is_better(f1, f2, e.fitness_scheme)
    return candidate1, f1
  else 
    return candidate2, f2
  end
end

function rank_by_fitness(e::Evaluator, candidates)
  # Note that we re-evaluate all candidates here. This might be wasteful and
  # we should cache if evaluations are costly.
  fitness = [(c[1], c[2], evaluate(e, c[1])) for c=candidates]
  sort(fitness; by = (t) -> t[3])
end

function fitness_is_within_ftol(e::Evaluator, ftolerance; index = 1)
  try
    fitness_is_within_ftol(e.problem, ftolerance, best_fitness(e.archive); index = index)
  catch
    false # In case we have no information yet about the best fitness
  end
end
