using BlackBoxOptim
using Distributions

abstract NaturalEvolutionStrategyOpt <: PopulationOptimizer

type SeparableNESOpt <: NaturalEvolutionStrategyOpt
  name::ASCIIString

  d::Integer                      # Number of dimensions
  lambda::Integer                 # Number of samples to take per iteration
  mu::Array{Float64,2}            # Average position for sampling each position
  sigma::Array{Float64,2}         # Average std deviation for sampling in position
  distr::Distribution             # Distribution we sample the step sizes from
  mu_learnrate::Float64
  sigma_learnrate::Float64
  last_s::Array{Float64,2}        # The s values sampled in the last call to ask
  population::Array{Float64,2}    # The last sampled values, now being evaluated

  SeparableNESOpt(searchSpace; lambda = false, mu_learnrate = 1.0, 
    sigma_learnrate = false) = begin

    numDimensions = length(mins(searchSpace))

    mu = rand(numDimensions, 1)
    sigma = ones(numDimensions, 1)
    distr = Normal(0, 1)

    lambda = lambda || convert(Int64, 4 + ceil(log(3*numDimensions)))
    sigma_learnrate = sigma_learnrate || calc_sigma_learnrate_for_snes(numDimensions)

    new("SeparableNES", numDimensions, lambda, mu, sigma, distr, 
      mu_learnrate, sigma_learnrate,
      eye(numDimensions), eye(numDimensions))

  end
end

NES_DefaultOptions = {
  "lambda" => false,          # If false it will be set based on the number of dimensions
  "mu_learnrate" => 1.0,
  "sigma_learnrate" => false, # If false it will be set based on the number of dimensions
}

function separable_nes(searchSpace; options = NES_DefaultOptions, population = false)
  # Note that we do not care about the population given as a parameter;
  # a NES generates a new population in every iteration so not used...
  SeparableNESOpt(searchSpace; lambda = options["lambda"], 
    mu_learnrate = options["mu_learnrate"], 
    sigma_learnrate = options["sigma_learnrate"])
end

calc_sigma_learnrate_for_snes(d) = (3 + log(d)) / (5 * sqrt(d))
calc_sigma_learnrate_for_nes(d) = (9 + 3 * log(d)) / (5 * d * sqrt(d))

# Get a set of new individuals to be ranked based on fitness.
function ask(snes::SeparableNESOpt)
  # Sample from N(0, 1)
  snes.last_s = s = randn(snes.d, snes.lambda)

  # Add in the mu and sigma's...
  #sampled_solutions = repmat(snes.mu, 1, l) + repmat(snes.sigma, 1, l) .* s
  # Quicker version as dimensions increase:
  snes.population = sampled_solutions = broadcast(+, snes.mu, broadcast(*, snes.sigma, s))

  # TODO!!! We must ensure all candidates are within the bounding box of the search space.

  # The rest of BlackBoxOptim still uses row-major order so transpose the 
  # individuals before returning them.
  mix_with_indices( sampled_solutions', 1:snes.lambda )
end

function mix_with_indices(candidates, indices = false)
  ary = Any[]
  for i in indices
    push!(ary, (candidates[i,:], i))
  end
  ary
end

# Tell the sNES the ranking of a set of candidates.
function tell!(snes::SeparableNESOpt, rankedCandidates)
  u = calc_utilities(rankedCandidates)'

  # Calc gradient
  gradient_mu = snes.last_s * u
  sq_s_minus1 = snes.last_s.^2 - 1
  gradient_sigma = sq_s_minus1 * u

  # Update the mean and sigma vectors based on the gradient
  snes.mu = snes.mu + snes.mu_learnrate * (snes.sigma .* gradient_mu)
  snes.sigma = snes.sigma .* exp(snes.sigma_learnrate / 2 * gradient_sigma)
end

function calc_utilities(rankedCandidates)
  num_candidates = length(rankedCandidates)

  # Second half has zero utility.
  treshold = convert(Int64, floor(num_candidates/2))
  second_half = zeros(num_candidates - treshold, 1)

  # While first half's utility decreases in linear steps 
  step_size = 1 / treshold
  first_half = linspace(1.0, step_size, treshold)

  # But the utilities should sum to 1 so we normalize, then return
  u = vcat(first_half, second_half)
  u = u ./ sum(u)

  # We must reorder the samples according to the order of the fitness in
  # rankedCandidates!!! Or we must reorder the utilities accordingly. The latter
  # is the preferred method and we can use the indices in rankedCandidates to 
  # accomplish it.
  u_ordered = zeros(num_candidates)
  for(i in 1:num_candidates)
    u_ordered[rankedCandidates[i][2]] = u[i]
  end

  u_ordered
end

#using BlackBoxOptim.Problems
#problem = "Sphere"
#p = BlackBoxOptim.Problems.examples[problem]
#fitness_for_opt(p, 2, 10, 10, separable_nes)