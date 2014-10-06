function test_pareto_dominance_func(pareto_dominance_func)
  test("simple cases with float vector fitness values of length 1") do
    @check pareto_dominance_func([0.0], [1.0])
    @check pareto_dominance_func([-100.0], [100.0])
    @check !pareto_dominance_func([-2.5], [-2.5])
    @check !pareto_dominance_func([12.32], [12.32])
  end

  test("simple cases with int vector fitness values of length 2") do
    @check pareto_dominance_func([0, 1], [1, 2])
    @check pareto_dominance_func([-10, -5], [0, 0])
    @check pareto_dominance_func([0, 1], [0, 2])
    @check pareto_dominance_func([-1, 2], [0, 2])

    @check !pareto_dominance_func([0, 1], [0, 1])
  end
end

describe("pareto_dominance") do
  describe("pareto_dominates_clear") do
    test_pareto_dominance_func(pareto_dominates_clear)
  end

  describe("pareto_dominates_fast") do
    test_pareto_dominance_func(pareto_dominates_fast)
  end
end