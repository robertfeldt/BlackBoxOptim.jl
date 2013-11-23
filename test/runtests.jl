include("helper.jl")

my_tests = [

  "utilities/test_latin_hypercube_sampling.jl",

  "test_fitness.jl",
  "test_population.jl",
  "test_bimodal_cauchy_distribution.jl",
  "test_search_space.jl",
  "test_frequency_adaptation.jl",

#  "test_random_search.jl",
  "test_differential_evolution.jl",
  "test_adaptive_differential_evolution.jl",
  "test_natural_evolution_strategies.jl",

  "problems/test_single_objective.jl",
]

for t in my_tests
  include(t)
end
