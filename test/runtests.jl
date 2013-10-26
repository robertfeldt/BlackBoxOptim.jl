using GlobalOptim
using FactCheck

my_tests = [
  "problems/single_objective.jl",
  "differential_evolution.jl"
]

for my_test in my_tests
  include(my_test)
end
