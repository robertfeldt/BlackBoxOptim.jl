using GlobalOptim
using FactCheck

my_tests = [
  "problems/single_objective.jl",
  "differential_evolution.jl"
]

println("Running tests:")

for my_test in my_tests
  #println(" * $(my_test)")
  FactCheck.TestSuite(my_test, my_test)
  include(my_test)
end
