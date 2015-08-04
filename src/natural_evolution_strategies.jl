using BlackBoxOptim
using Distributions

abstract NaturalEvolutionStrategyOpt <: PopulationOptimizer

type SeparableNESOpt <: NaturalEvolutionStrategyOpt
  d::Int                          # Number of dimensions
  lambda::Int                     # Number of samples to take per iteration
  mu::Vector{Float64}             # Average position for sampling each position
  sigma::Vector{Float64}          # Average std deviation for sampling in position
  distr::Distribution             # Distribution we sample the step sizes from
  mu_learnrate::Float64
  sigma_learnrate::Float64
  last_s::Array{Float64,2}        # The s values sampled in the last call to ask
  population::PopulationMatrix    # The last sampled values, now being evaluated
  utilities::Vector{Float64}      # The fitness shaping utility vector

  function SeparableNESOpt(searchSpace; lambda::Int = 0, mu_learnrate::Float64 = 1.0,
    sigma_learnrate::Float64 = 0.0)

    numDimensions = numdims(searchSpace)

    mu = rand(numDimensions)
    sigma = ones(numDimensions)
    distr = Normal(0, 1)

    if lambda == 0
      lambda = 4 + ceil(Int, log(3*numDimensions)) # default lambda
    end
    if sigma_learnrate == 0.0
      sigma_learnrate = calc_sigma_learnrate_for_snes(numDimensions) # default sigma learn rate
    end

    new(numDimensions, lambda, mu, sigma, distr,
      mu_learnrate, sigma_learnrate,
      zeros(numDimensions, lambda), zeros(numDimensions, lambda),
      # Most modern NES papers use log rather than linear fitness shaping.
      fitness_shaping_utilities_log(lambda))
  end
end


# We use a different ordering of the dimensions than other optimizers, so transpose.
population(o::NaturalEvolutionStrategyOpt) = o.population

const NES_DefaultOptions = @compat Dict{Symbol,Any}(
  :lambda => 0,              # If 0.0 it will be set based on the number of dimensions
  :mu_learnrate => 1.0,
  :sigma_learnrate => 0.0,   # If 0.0 it will be set based on the number of dimensions
)

function separable_nes(problem::OptimizationProblem, parameters)
  params = chain(NES_DefaultOptions, parameters)
  SeparableNESOpt(search_space(problem),
    lambda = params[:lambda],
    mu_learnrate = params[:mu_learnrate],
    sigma_learnrate = params[:sigma_learnrate])
end

calc_sigma_learnrate_for_snes(d) = (3 + log(d)) / (5 * sqrt(d))
calc_sigma_learnrate_for_nes(d) = (9 + 3 * log(d)) / (5 * d * sqrt(d))

# Get a set of new individuals to be ranked based on fitness.
function ask(snes::SeparableNESOpt)
  # Sample from N(0, 1)
  randn!(snes.last_s)

  # Add in the mu and sigma's...
  #sampled_solutions = repmat(snes.mu, 1, l) + repmat(snes.sigma, 1, l) .* s
  # Quicker version as dimensions increase:
  sampled_solutions = broadcast!(+, snes.population, snes.mu, broadcast(*, snes.sigma, snes.last_s))

  mix_with_indices(sampled_solutions, 1:snes.lambda)
end

function mix_with_indices(individuals::Matrix{Float64}, indices::Range)
  if popsize(individuals) != length(indices)
    throw(DimensionMismatch("The number of candidates does not match the number of indices"))
  end
  Candidate{Float64}[Candidate{Float64}(individuals[:,i], i) for i in indices]
end

# Tell the sNES the ranking of a set of candidates.
function tell!{F}(snes::SeparableNESOpt, rankedCandidates::Vector{Candidate{F}})
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

function assign_weights{F}(candidates::Vector{Candidate{F}}, u::Vector{Float64})
  n = length(candidates)
  # We must reorder the samples according to the order of the fitness in
  # candidates!!! Or we must reorder the utilities accordingly. The latter
  # is the preferred method and we can use the indices in candidates to
  # accomplish it.
  u[sortperm(candidates, by = x->x.fitness)]
end

# xNES is nice but scales badly with increasing dimension.
type XNESOpt <: NaturalEvolutionStrategyOpt
  d::Int                          # Number of dimensions
  lambda::Int                     # Number of samples to take per iteration
  utilities::Vector{Float64}      # Fitness utility to give to each rank
  x_learnrate::Float64
  a_learnrate::Float64
  A::Array{Float64,2}
  expA::Array{Float64,2}
  population::PopulationMatrix    # The last sampled values, now being evaluated
  x::Individual                   # The current incumbent (aka most likely value, mu etc)
  Z::Array{Float64,2}

  XNESOpt(searchSpace; lambda::Int = 0) = begin
    d = numdims(searchSpace)
    if lambda == 0
      lambda = 4 + 3*floor(Int, log(d))
    end
    x_learnrate = 1
    a_learnrate = 0.5 * minimum([1.0 / d, 0.25])
    x = rand_individual(searchSpace)

    new(d, lambda, fitness_shaping_utilities_log(lambda),
      x_learnrate, a_learnrate,
      zeros(d, d), zeros(d, d), zeros(d, lambda), x, zeros(d, lambda))
  end
end

function xnes(problem::OptimizationProblem, parameters)
  params = chain(NES_DefaultOptions, parameters)
  XNESOpt(search_space(problem); lambda = params[:lambda])
end

function ask(xnes::XNESOpt)
  xnes.expA = expm(xnes.A) # FIXME matrix allocation is expensive, use expm!() (when available)
  randn!(xnes.Z)
  broadcast!(+, xnes.population, xnes.x, (xnes.expA * xnes.Z))

  mix_with_indices(xnes.population, 1:xnes.lambda)
end

function tell!{F}(xnes::XNESOpt, rankedCandidates::Vector{Candidate{F}})
  u = assign_weights(rankedCandidates, xnes.utilities)

  # fixme improve memory footprint by using A_mul_B!() etc
  G = broadcast(*, u', xnes.Z) * xnes.Z' - eye(xnes.d)
  dx = xnes.x_learnrate * xnes.expA * (xnes.Z*u)
  dA = xnes.a_learnrate * G

  xnes.x += dx
  xnes.A += dA

  0
end

# Calculate the fitness shaping utilities vector using the log method.
function fitness_shaping_utilities_log(n::Int)
  u = [max(0.0, log(n / 2 + 1.0) - log(i)) for i in 1:n]
  u / sum(u)
end

# Calculate the fitness shaping utilities vector using the steps method.
function fitness_shaping_utilities_linear(n::Int)
  # Second half has zero utility.
  treshold = floor(Int, n/2)
  second_half = zeros(n - treshold)

  # While first half's utility decreases in linear steps
  step_size = 1 / treshold
  first_half = linspace(1.0, step_size, treshold)

  # But the utilities should sum to 1 so we normalize, then return
  u = vcat(first_half, second_half)
  u / sum(u)
end
