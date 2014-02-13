include("helper.jl")

my_slow_tests = [
  "problems/test_optimize_single_objective_problems.jl",
  "test_bboptimize.jl"
]

for t in my_slow_tests
  include(t)
end
