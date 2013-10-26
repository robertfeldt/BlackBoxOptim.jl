using GlobalOptim
using FactCheck

my_tests = [
  "problems/single_objective.jl",
  "differential_evolution.jl"
]

facts() do
  for my_test in my_tests
    FactCheck.TestSuite(my_test, my_test)
    include(my_test)
  end
end
