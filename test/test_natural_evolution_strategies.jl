using FactCheck

facts("sNES") do

function assign_weights_wrapper(candi_ixs::Vector{Int})
  candidates = BlackBoxOptim.Candidate{Float64}[BlackBoxOptim.Candidate{Float64}([0.0], i, NaN) for i in candi_ixs]
  u = BlackBoxOptim.fitness_shaping_utilities_linear(length(candi_ixs))
  BlackBoxOptim.assign_weights!(similar(u), candidates, u)
end

context("assign_weights") do
  context("when indices are already ordered") do
    u = assign_weights_wrapper([1, 2])
    @fact length(u) => 2
    @fact u[1] => 1.0
    @fact u[2] => 0.0

    u = assign_weights_wrapper([1, 2, 3])
    @fact length(u) => 3
    @fact u[1] => 1.0
    @fact u[2] => 0.0
    @fact u[3] => 0.0

    u = assign_weights_wrapper([1, 2, 3, 4])
    @fact length(u) => 4
    @fact u[1] => roughly(2/3)
    @fact u[2] => roughly(1/3)
    @fact u[3] => 0.0
    @fact u[4] => 0.0

    u = assign_weights_wrapper([1, 2, 3, 4, 5])
    @fact length(u) => 5
    @fact u[1] => roughly(2/3)
    @fact u[2] => roughly(1/3)
    @fact u[3] => 0.0
    @fact u[4] => 0.0
    @fact u[5] => 0.0

    u = assign_weights_wrapper([1, 2, 3, 4, 5, 6])
    @fact length(u) => 6
    @fact u[1] => roughly((3/3)/(3/3+2/3+1/3))
    @fact u[2] => roughly((2/3)/(3/3+2/3+1/3))
    @fact u[3] => roughly((1/3)/(3/3+2/3+1/3))
    @fact u[4] => 0.0
    @fact u[5] => 0.0
    @fact u[6] => 0.0
  end

  context("when indices are not ordered") do
    u = assign_weights_wrapper([4, 1, 2, 3])
    @fact u[4] => roughly(2/3)
    @fact u[1] => roughly(1/3)
    @fact u[2] => 0.0
    @fact u[3] => 0.0
  end
end

end
