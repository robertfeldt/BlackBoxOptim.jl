using GlobalOptim
using FactCheck

my_tests = [
  "problems/test_single_objective.jl",
  "test_differential_evolution.jl"
]

for my_test in my_tests
  include(my_test)
end
