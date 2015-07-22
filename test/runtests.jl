include("helper.jl")

my_tests = [

  "utilities/test_latin_hypercube_sampling.jl",
  "utilities/test_assign_ranks.jl",

  "test_parameters.jl",
  "test_fitness.jl",
  "test_population.jl",
  "test_bimodal_cauchy_distribution.jl",
  "test_search_space.jl",
  "test_mutation_operators.jl",
  "test_frequency_adaptation.jl",
  "test_archive.jl",

  "test_random_search.jl",
  "test_differential_evolution.jl",
  "test_adaptive_differential_evolution.jl",
  "test_natural_evolution_strategies.jl",

  "test_toplevel_bboptimize.jl",
  "test_smoketest_bboptimize.jl",

  "problems/test_single_objective.jl",

  "test_generating_set_search.jl",
  "test_direct_search_with_probabilistic_descent.jl",
]

for t in my_tests
  include(t)
end

exitstatus()
