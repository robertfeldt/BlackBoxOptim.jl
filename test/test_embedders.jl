facts("Embedding operators") do

context("RandomBound") do
  context("does nothing if within bounds") do
    @fact BlackBoxOptim.apply!(BlackBoxOptim.RandomBound([(0.0, 1.0)]),
                               [0.0], [0.0]) --> [0.0]

    @fact BlackBoxOptim.apply!(BlackBoxOptim.RandomBound([(0.0, 1.0), (10.0, 15.0)]),
                               [0.0, 11.4], [0.1, 12.3] ) --> [0.0, 11.4]
  end

  context("bounds if lower than min bound") do
    @fact BlackBoxOptim.apply!(BlackBoxOptim.RandomBound([(0.0, 1.0)]),
                               [-0.1], [0.0]) --> [0.0]

    res = BlackBoxOptim.apply!(BlackBoxOptim.RandomBound([(0.0, 1.0)]),
                               [-0.1], [0.5])
    @fact res[1] --> greater_than_or_equal(0.0)
    @fact res[1] --> less_than_or_equal(0.5)

    res = BlackBoxOptim.apply!(BlackBoxOptim.RandomBound([(-11.0, 1.0), (0.0, 1.0)]),
                               [-11.1, 0.5], [-10.8, 0.5])
    @fact (-10.8 >= res[1] >= -11.0) --> true
    @fact res[2] --> 0.5

    res = BlackBoxOptim.apply!(BlackBoxOptim.RandomBound([(30.0, 60.0), (-102.0, -1.0)]),
                               [50.4, -103.1], [49.6, -101.4])
    @fact res[1] --> 50.4
    @fact (-101.4 >= res[2] >= -102.0) --> true
  end

  context("bounds if higher than max bound") do
    @fact BlackBoxOptim.apply!(BlackBoxOptim.RandomBound([(0.0, 1.0)]),
                               [1.1], [1.0]) --> [1.0]

    res = BlackBoxOptim.apply!(BlackBoxOptim.RandomBound([(-10.0, 96.0)]),
                               [97.0], [95.0])
    @fact (95.0 <= res[1] <= 96.0) --> true
  end
end

end
