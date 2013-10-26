using GlobalOptim
using FactCheck

my_tests = [
  "problems/test_single_objective.jl",
  "test_differential_evolution.jl"
]

for t in my_tests
  include(t)
end
