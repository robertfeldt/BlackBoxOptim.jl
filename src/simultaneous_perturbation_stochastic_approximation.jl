abstract StochasticApproximationOptimizer <: Optimizer

SPSADefaultParameters = @compat Dict{Symbol,Any}(
  :Alpha => 0.602,  # The optimal value is 1.0 but values down to 0.602 often can give faster convergence
  :Gamma => 0.101,  # The optimal value is 1/6 but values down to 0.101 often can give faster convergence
  :a     => 0.0017,
  :c     => 1.9, # Recommendation is value 1.9 but that assumes noisy function, otherwise should be low
  :A     => 10
)

type SimultaneousPerturbationSA2{E<:EmbeddingOperator} <: StochasticApproximationOptimizer
  embed::E # embed candidate into search space
  parameters::Parameters
  k::Int64
  n::Int64
  theta::Individual
  delta_ck::Individual
end

function SimultaneousPerturbationSA2{E<:EmbeddingOperator}( embed::E, parameters )
    ss = parameters[:SearchSpace]
    n = numdims(ss)
    SimultaneousPerturbationSA2{E}(embed, Parameters(parameters, SPSADefaultParameters),
                                   0, n, rand_individual(ss), zeros(Float64, n))
end

# by default use RandomBound embedder
SimultaneousPerturbationSA2(parameters) = SimultaneousPerturbationSA2( RandomBound(parameters[:SearchSpace]), parameters )

name(spsa::SimultaneousPerturbationSA2) = "SPSA2 (Simultaneous Perturbation Stochastic Approximation, 1st order, 2 samples)"

sample_bernoulli_vector(n::Int) = 2.0 * round(rand(n)) - 1.0

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
  apply!(spsa.embed, theta_new, spsa.theta)
  spsa.theta = theta_new
  spsa.k += 1

  return 0
end
