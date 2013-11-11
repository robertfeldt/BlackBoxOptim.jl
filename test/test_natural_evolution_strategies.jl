using BlackBoxOptim.Problems
using FactCheck

facts("sNES") do

context("calc_utilities") do

  context("when indices are already ordered") do
    u = BlackBoxOptim.calc_utilities([(:dummy, 1.0), (:dummy, 2.0)])
    @fact length(u) => 2
    @fact u[1] => 1.0
    @fact u[2] => 0.0

    u = BlackBoxOptim.calc_utilities([(:d, 1.0), (:d, 2.0), (:d, 3.0)])
    @fact length(u) => 3
    @fact u[1] => 1.0
    @fact u[2] => 0.0
    @fact u[3] => 0.0

    u = BlackBoxOptim.calc_utilities([(:d, 1), (:d, 2), (:d, 3), (:d, 4)])
    @fact length(u) => 4
    @fact isapprox(u[1], 2/3) => true
    @fact isapprox(u[2], 1/3) => true
    @fact u[3] => 0.0
    @fact u[4] => 0.0

    u = BlackBoxOptim.calc_utilities([(:d, 1), (:d, 2), (:d, 3), (:d, 4), (:d, 5)])
    @fact length(u) => 5
    @fact isapprox(u[1], 2/3) => true
    @fact isapprox(u[2], 1/3) => true
    @fact u[3] => 0.0
    @fact u[4] => 0.0
    @fact u[5] => 0.0

    u = BlackBoxOptim.calc_utilities([(:d, 1), (:d, 2), (:d, 3), (:d, 4), (:d, 5), (:d, 6)])
    @fact length(u) => 6
    @fact isapprox(u[1], (((3/3)/(3/3+2/3+1/3)))) => true
    @fact isapprox(u[2], (((2/3)/(3/3+2/3+1/3)))) => true
    @fact isapprox(u[3], (((1/3)/(3/3+2/3+1/3)))) => true
    @fact u[4] => 0.0
    @fact u[5] => 0.0
    @fact u[6] => 0.0
  end

  context("when indices are not ordered") do
    u = BlackBoxOptim.calc_utilities([(:d, 4), (:d, 1), (:d, 2), (:d, 3)])
    @fact isapprox(u[4], 2/3) => true
    @fact isapprox(u[1], 1/3) => true
    @fact u[2] => 0.0
    @fact u[3] => 0.0
  end
end

#context("ask") do
#  p = BlackBoxOptim.Problems.examples["Sphere"]
#
#  problem = BlackBoxOptim.Problems.set_numdims!(2, p)
#
#  ss = search_space(problem)
#
#  opt = BlackBoxOptim.separable_nes(ss; population = false)
#
#  println("\n$(problem.name), n = $(numdims(problem)), optimizer = $(opt.name)")
#
#  best, fitness = BlackBoxOptim.run_optimizer_on_problem(opt, problem, 10)
#  fitness
#end

end