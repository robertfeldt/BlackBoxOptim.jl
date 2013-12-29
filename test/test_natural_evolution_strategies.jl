using FactCheck

facts("sNES") do

context("mix_with_indices") do
  cs = randn(2,3)
  r = BlackBoxOptim.mix_with_indices(cs)
  @fact length(r) => 2
  @fact r[1][2] => 1
  @fact r[2][2] => 2
  @fact r[1][1] => cs[1,:]
  @fact r[2][1] => cs[2,:]
end

assign_weights(rankedCandidates) = begin
  u = BlackBoxOptim.fitness_shaping_utilities_linear(length(rankedCandidates))
  BlackBoxOptim.assign_weights(rankedCandidates, u)
end

context("assign_weights") do

  context("when indices are already ordered") do
    u = assign_weights([(:dummy, 1.0), (:dummy, 2.0)])
    @fact length(u) => 2
    @fact u[1] => 1.0
    @fact u[2] => 0.0

    u = assign_weights([(:d, 1.0), (:d, 2.0), (:d, 3.0)])
    @fact length(u) => 3
    @fact u[1] => 1.0
    @fact u[2] => 0.0
    @fact u[3] => 0.0

    u = assign_weights([(:d, 1), (:d, 2), (:d, 3), (:d, 4)])
    @fact length(u) => 4
    @fact isapprox(u[1], 2/3) => true
    @fact isapprox(u[2], 1/3) => true
    @fact u[3] => 0.0
    @fact u[4] => 0.0

    u = assign_weights([(:d, 1), (:d, 2), (:d, 3), (:d, 4), (:d, 5)])
    @fact length(u) => 5
    @fact isapprox(u[1], 2/3) => true
    @fact isapprox(u[2], 1/3) => true
    @fact u[3] => 0.0
    @fact u[4] => 0.0
    @fact u[5] => 0.0

    u = assign_weights([(:d, 1), (:d, 2), (:d, 3), (:d, 4), (:d, 5), (:d, 6)])
    @fact length(u) => 6
    @fact isapprox(u[1], (((3/3)/(3/3+2/3+1/3)))) => true
    @fact isapprox(u[2], (((2/3)/(3/3+2/3+1/3)))) => true
    @fact isapprox(u[3], (((1/3)/(3/3+2/3+1/3)))) => true
    @fact u[4] => 0.0
    @fact u[5] => 0.0
    @fact u[6] => 0.0
  end

  context("when indices are not ordered") do
    u = assign_weights([(:d, 4), (:d, 1), (:d, 2), (:d, 3)])
    @fact isapprox(u[4], 2/3) => true
    @fact isapprox(u[1], 1/3) => true
    @fact u[2] => 0.0
    @fact u[3] => 0.0
  end
end

end