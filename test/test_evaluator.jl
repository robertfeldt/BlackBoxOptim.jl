function evaluator_tests(make_eval::Function)
  context("Basic evaluation with a single-objective function to be minimized") do
    e = make_eval()
    @fact BlackBoxOptim.num_evals(e) --> 0

    fit1 = fitness([0.0, 1.0], e)
    @fact BlackBoxOptim.num_evals(e) --> 1
    @fact BlackBoxOptim.last_fitness(e) --> fit1

    fit2 = fitness([2.0, 1.0], e)
    @fact BlackBoxOptim.num_evals(e) --> 2
    @fact BlackBoxOptim.last_fitness(e) --> fit2

    @fact BlackBoxOptim.numdims(e) --> 2

    @fact BlackBoxOptim.is_better(0.0, 1.0, e) --> true
    @fact BlackBoxOptim.is_better(0.0, -1.0, e) --> false

    a, b = [0.0, 0.5], [1.0, 1.0]

    @fact BlackBoxOptim.is_better(a, 1.0, e) --> true
    @fact BlackBoxOptim.last_fitness(e) --> 0.25

    @fact BlackBoxOptim.is_better(b, 1.0, e) --> false
    @fact BlackBoxOptim.last_fitness(e) --> 2.0

    @fact BlackBoxOptim.best_of(a, b, e) --> (a, 0.25)
    @fact BlackBoxOptim.best_of(b, a, e) --> (a, 0.25)
    BlackBoxOptim.shutdown!(e)
  end

  context("update_fitness!()") do
    e = make_eval()

    candidates = BlackBoxOptim.Candidate{Float64}[BlackBoxOptim.Candidate{Float64}(clamp(randn(2), -1.0, 1.0)) for i in 1:10]
    BlackBoxOptim.update_fitness!(e, candidates)
    @fact BlackBoxOptim.num_evals(e) --> length(candidates)
    @fact all(c -> isfinite(c.fitness), candidates) --> true
    BlackBoxOptim.shutdown!(e)
  end

  context("rank_by_fitness!()") do
    e = make_eval()

    candidates = BlackBoxOptim.Candidate{Float64}[BlackBoxOptim.Candidate{Float64}(clamp(randn(2), -1.0, 1.0)) for i in 1:10]
    # partially evaluate fitness
    BlackBoxOptim.update_fitness!(e, candidates[1:5])
    @fact BlackBoxOptim.num_evals(e) --> 5
    # complete fitness evaluation and sort by it
    BlackBoxOptim.rank_by_fitness!(e, candidates)
    @fact BlackBoxOptim.num_evals(e) --> 10
    @fact sortperm(candidates, by = fitness) --> collect(1:10)
    BlackBoxOptim.shutdown!(e)
  end
end

facts("Evaluator") do
  # Set up a small example problem
  f(x) = sum(x.^2)
  p = minimization_problem(f, "", (-1.0, 1.0), 2)
  context("ProblemEvaluator") do
    evaluator_tests(() -> BlackBoxOptim.ProblemEvaluator(p))
  end

if BlackBoxOptim.enable_parallel_methods
  context("ParallelEvaluator") do
    evaluator_tests(() -> BlackBoxOptim.ParallelEvaluator(p, pids=workers()))
  end
end

end
