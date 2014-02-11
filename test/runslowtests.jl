include("helper.jl")

my_slow_tests = [
  "problems/test_optimize_single_objective_problems.jl",
  "test_bboptimize.jl",
  "test_compare_optimizers.jl",
]

for t in my_slow_tests
  include(t)
end
