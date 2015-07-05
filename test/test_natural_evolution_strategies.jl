using FactCheck

facts("sNES") do

context("mix_with_indices") do
  cs = randn(2,3)
  @fact_throws DimensionMismatch BlackBoxOptim.mix_with_indices(cs, 1:2)
  r = BlackBoxOptim.mix_with_indices(cs, 1:3)
  @fact length(r) => 3
  @fact r[1].index => 1
  @fact r[2].index => 2
  @fact r[3].index => 3
  @fact r[1].params => cs[:,1]
  @fact r[2].params => cs[:,2]
  @fact r[3].params => cs[:,3]
end

function assign_weights{F}(fitnesses::Vector{F})
  candidates = [BlackBoxOptim.Candidate{F}([0.0], -1, f) for f in fitnesses]
  u = BlackBoxOptim.fitness_shaping_utilities_linear(length(candidates))
  BlackBoxOptim.assign_weights(candidates, u)
end

context("assign_weights") do
  context("when indices are already ordered") do
    u = assign_weights([1.0, 2.0])
    @fact length(u) => 2
    @fact u[1] => 1.0
    @fact u[2] => 0.0

    u = assign_weights([1.0, 2.0, 3.0])
    @fact length(u) => 3
    @fact u[1] => 1.0
    @fact u[2] => 0.0
    @fact u[3] => 0.0

    u = assign_weights([1, 2, 3, 4])
    @fact length(u) => 4
    @fact u[1] => roughly(2/3)
    @fact u[2] => roughly(1/3)
    @fact u[3] => 0.0
    @fact u[4] => 0.0

    u = assign_weights([1, 2, 3, 4, 5])
    @fact length(u) => 5
    @fact u[1] => roughly(2/3)
    @fact u[2] => roughly(1/3)
    @fact u[3] => 0.0
    @fact u[4] => 0.0
    @fact u[5] => 0.0

    u = assign_weights([1, 2, 3, 4, 5, 6])
    @fact length(u) => 6
    @fact u[1] => roughly((3/3)/(3/3+2/3+1/3))
    @fact u[2] => roughly((2/3)/(3/3+2/3+1/3))
    @fact u[3] => roughly((1/3)/(3/3+2/3+1/3))
    @fact u[4] => 0.0
    @fact u[5] => 0.0
    @fact u[6] => 0.0
  end

  context("when indices are not ordered") do
    u = assign_weights([4, 1, 2, 3])
    @fact u[4] => roughly(2/3)
    @fact u[1] => roughly(1/3)
    @fact u[2] => 0.0
    @fact u[3] => 0.0
  end
end

end
