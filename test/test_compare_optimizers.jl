facts("compare_optimizers") do
  context("comparing optimizers on more than one problem") do
    #methods = BlackBoxOptim.MethodNames
    methods = [:adaptive_de_rand_1_bin_radiuslimited, :generating_set_search, :separable_nes, :xnes
      #:random_search,
      #:simultaneous_perturbation_stochastic_approximation
    ]
    problems = BlackBoxOptim.as_fixed_dim_problem_set(BlackBoxOptim.example_problems, [50, 500])
    ranks, fitnesses = BlackBoxOptim.compare_optimizers(problems; max_time = 60.0, methods = methods)
    @fact size(ranks, 2) --> length(problems)
  end
end