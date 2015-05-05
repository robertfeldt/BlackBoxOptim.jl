using BlackBoxOptim
using Distributions

abstract NaturalEvolutionStrategyOpt <: PopulationOptimizer

type SeparableNESOpt <: NaturalEvolutionStrategyOpt
  d::Integer                      # Number of dimensions
  lambda::Integer                 # Number of samples to take per iteration
  mu::Array{Float64,2}            # Average position for sampling each position
  sigma::Array{Float64,2}         # Average std deviation for sampling in position
  distr::Distribution             # Distribution we sample the step sizes from
  mu_learnrate::Float64
  sigma_learnrate::Float64
  last_s::Array{Float64,2}        # The s values sampled in the last call to ask
  population::Array{Float64,2}    # The last sampled values, now being evaluated
  utilities::Array{Float64,2}     # The fitness shaping utility vector

  SeparableNESOpt(searchSpace; lambda = false, mu_learnrate = 1.0, 
    sigma_learnrate = false) = begin

    numDimensions = length(mins(searchSpace))

    mu = rand(numDimensions, 1)
    sigma = ones(numDimensions, 1)
    distr = Normal(0, 1)

    lambda = lambda || convert(Int, 4 + ceil(log(3*numDimensions)))
    sigma_learnrate = sigma_learnrate || calc_sigma_learnrate_for_snes(numDimensions)

    new(numDimensions, lambda, mu, sigma, distr, 
      mu_learnrate, sigma_learnrate,
      eye(numDimensions), eye(numDimensions), 
      # Most modern NES papers use log rather than linear fitness shaping.
      fitness_shaping_utilities_log(lambda))
  end
end


# We use a different ordering of the dimensions than other optimizers, so transpose.
population(o::NaturalEvolutionStrategyOpt) = o.population'

NES_DefaultOptions = {
  "lambda" => false,          # If false it will be set based on the number of dimensions
  "mu_learnrate" => 1.0,
  "sigma_learnrate" => false, # If false it will be set based on the number of dimensions
}

function separable_nes(parameters)
  params = mergeparam(NES_DefaultOptions, parameters)
  SeparableNESOpt(params[:SearchSpace]; 
    lambda = params["lambda"], 
    mu_learnrate = params["mu_learnrate"], 
    sigma_learnrate = params["sigma_learnrate"])
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

  # The rest of BlackBoxOptim still uses row-major order so transpose the 
  # individuals before returning them.
  mix_with_indices( sampled_solutions', 1:snes.lambda )
end

function mix_with_indices(candidates, indices = false)
  if indices == false
    indices = 1:size(candidates,1)
  end
  ary = Any[]
  for i in indices
    push!(ary, (candidates[i,:], i))
  end
  ary
end

# Tell the sNES the ranking of a set of candidates.
function tell!(snes::SeparableNESOpt, rankedCandidates)
  u = assign_weights(rankedCandidates, snes.utilities)

  # Calc gradient
  gradient_mu = snes.last_s * u
  sq_s_minus1 = snes.last_s.^2 - 1
  gradient_sigma = sq_s_minus1 * u

  # Update the mean and sigma vectors based on the gradient
  snes.mu = snes.mu + snes.mu_learnrate * (snes.sigma .* gradient_mu)
  snes.sigma = snes.sigma .* exp(snes.sigma_learnrate / 2 * gradient_sigma)

  # There is no notion of how many was better in NES so return 0
  0
end

function assign_weights(rankedCandidates, u)
  n = length(rankedCandidates)
  # We must reorder the samples according to the order of the fitness in
  # rankedCandidates!!! Or we must reorder the utilities accordingly. The latter
  # is the preferred method and we can use the indices in rankedCandidates to 
  # accomplish it.
  u_ordered = zeros(n, 1)
  for(i in 1:n)
    u_ordered[rankedCandidates[i][2]] = u[i]
  end

  u_ordered
end

# xNES is nice but scales badly with increasing dimension.
type XNESOpt <: NaturalEvolutionStrategyOpt
  d::Integer                      # Number of dimensions
  lambda::Integer                 # Number of samples to take per iteration
  utilities::Array{Float64,2}     # Fitness utility to give to each rank
  x_learnrate::Float64
  a_learnrate::Float64
  A::Array{Float64,2}
  expA::Array{Float64,2}
  population::Array{Float64,2}    # The last sampled values, now being evaluated
  x::Array{Float64,2}             # The current incumbent (aka most likely value, mu etc)
  Z::Array{Float64,2}

  XNESOpt(searchSpace; lambda = false) = begin
    d = numdims(searchSpace)
    lambda = lambda || convert(Int, 4 + 3*floor(log(d)))
    x_learnrate = 1
    a_learnrate = 0.5 * minimum([1.0 / d, 0.25])
    x = rand_individual(searchSpace)

    new(d, lambda, fitness_shaping_utilities_log(lambda), 
      x_learnrate, a_learnrate, 
      zeros(d, d), zeros(d, d), zeros(d, lambda), x, zeros(d, d))
  end
end

function xnes(parameters)
  params = mergeparam(NES_DefaultOptions, parameters)
  XNESOpt(params[:SearchSpace]; lambda = params["lambda"])
end

function ask(xnes::XNESOpt)
  xnes.expA = expm(xnes.A)
  xnes.Z = randn(xnes.d, xnes.lambda)
  xnes.population = broadcast(+, xnes.x, (xnes.expA * xnes.Z))

  # The rest of BlackBoxOptim still uses row-major order so transpose the 
  # individuals before returning them.
  mix_with_indices( xnes.population', 1:xnes.lambda )
end

function tell!(xnes::XNESOpt, rankedCandidates)
  u = assign_weights(rankedCandidates, xnes.utilities)

  G = (repmat(u', xnes.d, 1) .* xnes.Z) * xnes.Z' - eye(xnes.d)
  dx = xnes.x_learnrate * xnes.expA * (xnes.Z*u)
  dA = xnes.a_learnrate * G
  
  xnes.x += dx
  xnes.A += dA

  0
end

# Calculate the fitness shaping utilities vector using the log method.
function fitness_shaping_utilities_log(n)
  u = maximum(hcat(zeros(n, 1), log(n / 2 + 1.0) - log(1:n)), 2)
  u / sum(u)
end

# Calculate the fitness shaping utilities vector using the steps method.
function fitness_shaping_utilities_linear(n)
  # Second half has zero utility.
  treshold = convert(Int, floor(n/2))
  second_half = zeros(n - treshold, 1)

  # While first half's utility decreases in linear steps 
  step_size = 1 / treshold
  first_half = linspace(1.0, step_size, treshold)

  # But the utilities should sum to 1 so we normalize, then return
  u = vcat(first_half, second_half)
  u / sum(u)  
end
