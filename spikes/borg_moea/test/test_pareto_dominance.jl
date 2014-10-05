describe("pareto_dominance") do
  test("simple cases with float vector fitness values of length 1") do
    @check pareto_dominates([0.0], [1.0])
  end

  test("simple cases with int vector fitness values of length 2") do
    @check pareto_dominates([0, 1], [1, 2])
    @check pareto_dominates([0, 1], [0, 2])
    @check !pareto_dominates([0, 1], [0, 1])
  end
end