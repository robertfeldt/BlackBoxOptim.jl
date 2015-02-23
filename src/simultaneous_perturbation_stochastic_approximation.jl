abstract StochasticApproximationOptimizer <: Optimizer

SPSADefaultParameters = {
  :Alpha => 0.602,  # The optimal value is 1.0 but values down to 0.602 often can give faster convergence
  :Gamma => 0.101,  # The optimal value is 1/6 but values down to 0.101 often can give faster convergence
  :a     => 0.0017,
  :c     => 1.9, # Recommendation is value 1.9 but that assumes noisy function, otherwise should be low
  :A     => 10
}

type SimultaneousPerturbationSA2 <: StochasticApproximationOptimizer
  search_space::SearchSpace
  parameters::Parameters
  k::Int
  n::Int
  theta::Array{Float64, 2}
  delta_ck::Array{Float64, 2}

  SimultaneousPerturbationSA2(parameters) = begin
    ss = parameters[:SearchSpace]
    n = numdims(ss)
    new(ss, Parameters(parameters, SPSADefaultParameters),
      0, n, rand_individual(ss), zeros(Float64, n, 2))
  end
end

name(spsa::SimultaneousPerturbationSA2) = "SPSA2 (Simultaneous Perturbation Stochastic Approximation, 1st order, 2 samples)"

sample_bernoulli_vector(n::Int) = 2.0 * round(rand(n,1)) - 1.0

function ask(spsa::SimultaneousPerturbationSA2)

  delta = sample_bernoulli_vector(spsa.n)
  ck = spsa.parameters[:c]/(spsa.k + 1)^spsa.parameters[:Gamma]
  spsa.delta_ck = ck * delta

  theta_plus = spsa.theta + spsa.delta_ck
  theta_minus = spsa.theta - spsa.delta_ck

  [(theta_plus, 1), (theta_minus, 2)]

end

function tell!(spsa::SimultaneousPerturbationSA2, rankedCandidates)

  # Use index of rank to get right values for yplus and yminus, respectively.
  if rankedCandidates[1][2] == 1
    yplus = rankedCandidates[1][3]
    yminus = rankedCandidates[2][3]
  else
    yplus = rankedCandidates[2][3]
    yminus = rankedCandidates[1][3]
  end

  # Estimate gradient.
  ghat = (yplus - yminus) ./ (2.0 * spsa.delta_ck)

  # Calc new estimate of theta based on estimate of gradient, ghat.
  ak = spsa.parameters[:a]/(spsa.k + 1 + spsa.parameters[:A])^spsa.parameters[:Alpha]
  theta_new = spsa.theta - ak * ghat
  rand_bound_from_target!(theta_new, spsa.theta, spsa.search_space)
  spsa.theta = theta_new
  spsa.k += 1

  return 0
end
