using GlobalOptim.Problems

facts("Sphere") do
  p = GlobalOptim.Problems.examples["Sphere"]
  sphere = p.funcs[1]

  @fact sphere([0]) => 0

  @fact sphere([1]) => 1

  @fact sphere([1, 2]) => 5

  @fact sphere([1, 2, 3]) => 14

  @fact sphere([-1, 2, -3]) => 14

  @fact_throws sphere([])

  @fact numdims(p) => 2 # Default is 2 if not specified
  @fact search_space(p) => [(-100.0, 100.0), (-100.0, 100.0)]
end

facts("Schwefel2.22") do
  p = GlobalOptim.Problems.examples["Schwefel2.22"]
  schwefel2_22 = p.funcs[1]

  @fact schwefel2_22([0]) => 0

  @fact schwefel2_22([1]) => 2

  @fact schwefel2_22([1, 2]) => (1+2)+(1*2)

  @fact schwefel2_22([1, 2, 3]) => (1+2+3)+(1*2*3)

  @fact schwefel2_22([-1, 2, -3]) => (1+2+3)+(1*2*3)

  @fact_throws schwefel2_22([])

  @fact numdims(p) => 2 # Default is 2 if not specified
  @fact search_space(p) => [(-10.0, 10.0), (-10.0, 10.0)]
end

facts("Schwefel1.2") do
  p = GlobalOptim.Problems.examples["Schwefel1.2"]
  schwefel1_2 = p.funcs[1]

  @fact schwefel1_2([0]) => 0

  @fact schwefel1_2([1]) => 1

  @fact schwefel1_2([1, 2]) => 1+9

  @fact schwefel1_2([1, 2, 3]) => 1+9+36

  @fact schwefel1_2([-1, 2, -3]) => 1+1+4

  @fact schwefel1_2([]) => 0

  @fact numdims(p) => 2 # Default is 2 if not specified
  @fact search_space(p) => [(-100.0, 100.0), (-100.0, 100.0)]
end

facts("Schwefel2.21") do
  p = GlobalOptim.Problems.examples["Schwefel2.21"]
  schwefel2_21 = p.funcs[1]

  @fact schwefel2_21([0]) => 1

  @fact schwefel2_21([1]) => 1

  @fact schwefel2_21([1, 2]) => 2

  @fact schwefel2_21([1, 2, 3]) => 3

  @fact schwefel2_21([-1, 2, -3]) => 3

  @fact_throws schwefel2_21([])

  @fact numdims(p) => 2 # Default is 2 if not specified
  @fact search_space(p) => [(-100.0, 100.0), (-100.0, 100.0)]
end
