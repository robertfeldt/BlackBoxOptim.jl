"""
`AskTellOptimizer` that utilizes randomization to generate the candidates.
"""
abstract type StochasticApproximationOptimizer <: AskTellOptimizer end

const SPSADefaultParameters = ParamsDict(
    :Alpha => 0.602,  # The optimal value is 1.0 but values down to 0.602 often can give faster convergence
    :Gamma => 0.101,  # The optimal value is 1/6 but values down to 0.101 often can give faster convergence
    :a     => 0.0017,
    :c     => 1.9, # Recommendation is value 1.9 but that assumes noisy function, otherwise should be low
    :A     => 10
)

mutable struct SimultaneousPerturbationSA2{E<:EmbeddingOperator} <: StochasticApproximationOptimizer
    embed::E # embed candidate into search space
    parameters::Parameters
    k::Int
    n::Int
    theta::Individual
    delta_ck::Individual

    function SimultaneousPerturbationSA2(problem::OptimizationProblem, embed::E, parameters::Parameters) where {E<:EmbeddingOperator}
        ss = search_space(problem)
        n = numdims(ss)
        new{E}(embed, chain(SPSADefaultParameters, parameters),
               0, n, rand_individual(ss), zeros(Float64, n))
    end
end

# by default use RandomBound embedder
SimultaneousPerturbationSA2(problem::OptimizationProblem, parameters::Parameters) =
    SimultaneousPerturbationSA2(problem, RandomBound(search_space(problem)), parameters)

name(spsa::SimultaneousPerturbationSA2) = "SPSA2 (Simultaneous Perturbation Stochastic Approximation, 1st order, 2 samples)"

sample_bernoulli_vector(n::Int) = [2.0*round(rand()) - 1.0 for _ in 1:n]

function ask(spsa::SimultaneousPerturbationSA2)
    delta = sample_bernoulli_vector(spsa.n)
    ck = spsa.parameters[:c]/(spsa.k + 1)^spsa.parameters[:Gamma]
    spsa.delta_ck = ck * delta

    theta_plus = spsa.theta + spsa.delta_ck
    theta_minus = spsa.theta - spsa.delta_ck

    [Candidate{Float64}(theta_plus, 1),
     Candidate{Float64}(theta_minus, 2)]
end

function tell!(spsa::SimultaneousPerturbationSA2, rankedCandidates::AbstractVector{<:Candidate})
    # Use index of rank to get right values for yplus and yminus, respectively.
    if rankedCandidates[1].index == 1
        yplus = rankedCandidates[1].fitness
        yminus = rankedCandidates[2].fitness
    else
        yplus = rankedCandidates[2].fitness
        yminus = rankedCandidates[1].fitness
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
