using GlobalOptim
using FactCheck

my_tests = ["problems/single_objective.jl"]

println("Running tests:")

for my_test in my_tests
    println(" * $(my_test)")
    include(my_test)
end