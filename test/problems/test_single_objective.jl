facts("Single objective functions") do

context("Sphere") do
  p = BlackBoxOptim.example_problems["Sphere"]
  sphere = objfunc(p)

  @fact sphere([0]) --> 0

  @fact sphere([1]) --> 1

  @fact sphere([1, 2]) --> 5

  @fact sphere([1, 2, 3]) --> 14

  @fact sphere([-1, 2, -3]) --> 14

  @fact_throws sphere([])

  p2 = fixed_dim_problem(p, 3)
  @fact numdims(p2) --> 3
  @fact ranges(search_space(p2)) --> [(-100.0, 100.0), (-100.0, 100.0), (-100.0, 100.0)]
end

context("Schwefel2.22") do
  p = BlackBoxOptim.example_problems["Schwefel2.22"]
  schwefel2_22 = objfunc(p)

  @fact schwefel2_22([0]) --> 0

  @fact schwefel2_22([1]) --> 2

  @fact schwefel2_22([1, 2]) --> (1+2)+(1*2)

  @fact schwefel2_22([1, 2, 3]) --> (1+2+3)+(1*2*3)

  @fact schwefel2_22([-1, 2, -3]) --> (1+2+3)+(1*2*3)

  @fact_throws schwefel2_22([])

  p2 = fixed_dim_problem(p, 4)
  @fact numdims(p2) --> 4
  @fact ranges(search_space(p2)) --> [(-10.0, 10.0), (-10.0, 10.0), (-10.0, 10.0), (-10.0, 10.0)]
end

context("Schwefel1.2") do
  p = BlackBoxOptim.example_problems["Schwefel1.2"]
  schwefel1_2 = objfunc(p)

  @fact schwefel1_2([0]) --> 0

  @fact schwefel1_2([1]) --> 1

  @fact schwefel1_2([1, 2]) --> 1+9

  @fact schwefel1_2([1, 2, 3]) --> 1+9+36

  @fact schwefel1_2([-1, 2, -3]) --> 1+1+4

  @fact schwefel1_2([]) --> 0
end

context("Schwefel2.21") do
  p = BlackBoxOptim.example_problems["Schwefel2.21"]
  schwefel2_21 = objfunc(p)

  @fact schwefel2_21([0]) --> 0

  @fact schwefel2_21([1]) --> 1

  @fact schwefel2_21([1, 2]) --> 2

  @fact schwefel2_21([1, 2, 3]) --> 3

  @fact schwefel2_21([-1, 2, -3]) --> 3

  @fact_throws schwefel2_21([])
end

context("Rosenbrock") do
  p = BlackBoxOptim.example_problems["Rosenbrock"]
  rosenbrock = objfunc(p)

  @fact rosenbrock([1, 2]) --> 100

  @fact rosenbrock([1, 2, 3]) --> 201

  @fact rosenbrock([-1, 2, -3]) --> 5005

  @fact_throws rosenbrock([])
end

end