using FactCheck
using GlobalOptim.SingleObjectiveProblems

facts("Sphere") do
  sphere = GlobalOptim.SingleObjectiveProblems.examples["Sphere"].f

  @fact sphere([0]) => 0

  @fact sphere([1]) => 1

  @fact sphere([1, 2]) => 5

  @fact sphere([1, 2, 3]) => 14

  @fact sphere([-1, 2, -3]) => 14

  @fact_throws sphere([])
end

facts("Schwefel2_22") do
  schwefel2_22 = GlobalOptim.SingleObjectiveProblems.examples["Schwefel2.22"].f

  @fact schwefel2_22([0]) => 0

  @fact schwefel2_22([1]) => 2

  @fact schwefel2_22([1, 2]) => (1+2)+(1*2)

  @fact schwefel2_22([1, 2, 3]) => (1+2+3)+(1*2*3)

  @fact schwefel2_22([-1, 2, -3]) => (1+2+3)+(1*2*3)

  @fact_throws schwefel2_22([])
end

facts("Schwefel1_2") do
  schwefel1_2 = GlobalOptim.SingleObjectiveProblems.examples["Schwefel1.2"].f

  @fact schwefel1_2([0]) => 0

  @fact schwefel1_2([1]) => 1

  @fact schwefel1_2([1, 2]) => 1+9

  @fact schwefel1_2([1, 2, 3]) => 1+9+36

  @fact schwefel1_2([-1, 2, -3]) => 1+1+4

  @fact schwefel1_2([]) => 0
end

facts("Schwefel2_21") do
  schwefel2_21 = GlobalOptim.SingleObjectiveProblems.examples["Schwefel2.21"].f

  @fact schwefel2_21([0]) => 1

  @fact schwefel2_21([1]) => 1

  @fact schwefel2_21([1, 2]) => 2

  @fact schwefel2_21([1, 2, 3]) => 3

  @fact schwefel2_21([-1, 2, -3]) => 3

  @fact_throws schwefel2_21([])
end
