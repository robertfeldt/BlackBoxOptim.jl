facts("compare_optimizers") do
  context("comparing optimizers on more than one problem") do
    problems = BlackBoxOptim.as_fixed_dim_problem_set(BlackBoxOptim.example_problems, 30)
    ranks, fitnesses = BlackBoxOptim.compare_optimizers(problems; max_time = 5.0)
    @fact size(ranks, 1) => length(problems)
  end
end