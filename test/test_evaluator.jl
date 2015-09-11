facts("Evaluator") do

  # Set up a small example problem
  f(x) = sum(x.^2)
  P = BlackBoxOptim.fixeddim_problem(f; dims = 2)

  context("Basic evaluation with a single-objective function to be minimized") do

    e = ProblemEvaluator(P)

    @fact BlackBoxOptim.num_evals(e) --> 0

    f = BlackBoxOptim.evaluate(e, [0.0, 1.0])
    @fact BlackBoxOptim.num_evals(e) --> 1
    @fact BlackBoxOptim.last_fitness(e) --> f

    f2 = BlackBoxOptim.evaluate(e, [2.0, 1.0])
    @fact BlackBoxOptim.num_evals(e) --> 2
    @fact BlackBoxOptim.last_fitness(e) --> f2

    @fact BlackBoxOptim.numdims(e) --> 2

    @fact BlackBoxOptim.is_better(e, 0.0, 1.0) --> true
    @fact BlackBoxOptim.is_better(e, 0.0, -1.0) --> false

    a, b = [0, 0.5], [1, 1]

    @fact BlackBoxOptim.is_better(e, a, 1.0) --> true
    @fact BlackBoxOptim.last_fitness(e) --> 0.25

    @fact BlackBoxOptim.is_better(e, b, 1.0) --> false
    @fact BlackBoxOptim.last_fitness(e) --> 2.0

    @fact BlackBoxOptim.best_of(e, a, b) --> (a, 0.25)
    @fact BlackBoxOptim.best_of(e, b, a) --> (a, 0.25)

  end

end