abstract NaturalEvolutionStrategyOpt <: PopulationOptimizer

type SeparableNESOpt{F,E<:EmbeddingOperator} <: NaturalEvolutionStrategyOpt
  embed::E                        # operator embedding into the search space
  lambda::Int                     # Number of samples to take per iteration
  mu::Vector{Float64}             # Average position for sampling each position
  sigma::Vector{Float64}          # Average std deviation for sampling in position
  distr::Distribution             # Distribution we sample the step sizes from
  mu_learnrate::Float64
  sigma_learnrate::Float64
  last_s::Array{Float64,2}        # The s values sampled in the last call to ask
  candidates::Vector{Candidate{F}}# The last sampled values, now being evaluated
  utilities::Vector{Float64}      # The fitness shaping utility vector

  function SeparableNESOpt(embed::E;
    lambda::Int = 0,
    mu_learnrate::Float64 = 1.0,
    sigma_learnrate::Float64 = 0.0,
    ini_x = nothing)
    d = numdims(search_space(embed))

    if lambda == 0
      lambda = 4 + ceil(Int, log(3*d)) # default lambda
    end
    if sigma_learnrate == 0.0
      sigma_learnrate = calc_sigma_learnrate_for_snes(d) # default sigma learn rate
    end
    if ini_x === nothing
      ini_x = rand_individual(search_space(embed))
    else
      ini_x = copy(ini_x::Individual)
    end

    new(embed, lambda,
      ini_x, ones(d),
      Normal(0, 1),
      mu_learnrate, sigma_learnrate,
      zeros(d, lambda),
      Candidate{F}[Candidate{F}(Array(Float64, d), i) for i in 1:lambda],
      # Most modern NES papers use log rather than linear fitness shaping.
      fitness_shaping_utilities_log(lambda))
  end
end

# We use a different ordering of the dimensions than other optimizers, so transpose.
population(o::NaturalEvolutionStrategyOpt) = o.candidates

const NES_DefaultOptions = @compat Dict{Symbol,Any}(
  :lambda => 0,              # If 0 it will be set based on the number of dimensions
  :ini_x => nothing,         # starting point, "nothing" generates random point in a search space
  :mu_learnrate => 1.0,
  :sigma_learnrate => 0.0,   # If 0.0 it will be set based on the number of dimensions
)

function separable_nes(problem::OptimizationProblem, parameters)
  params = chain(NES_DefaultOptions, parameters)
  embed = RandomBound(search_space(problem))
  SeparableNESOpt{fitness_type(problem), typeof(embed)}(embed,
    lambda = params[:lambda],
    mu_learnrate = params[:mu_learnrate],
    sigma_learnrate = params[:sigma_learnrate],
    ini_x = params[:ini_x])
end

calc_sigma_learnrate_for_snes(d) = (3 + log(d)) / (5 * sqrt(d))
calc_sigma_learnrate_for_nes(d) = (9 + 3 * log(d)) / (5 * d * sqrt(d))

# Get a set of new individuals to be ranked based on fitness.
function ask(snes::SeparableNESOpt)
  # Sample from N(0, 1)
  randn!(snes.last_s)

  for i in eachindex(snes.candidates)
    candi = snes.candidates[i]
    @inbounds for j in 1:length(candi.params)
      candi.params[j] = snes.mu[j] + snes.sigma[j] * snes.last_s[j, i]
    end
    apply!(snes.embed, candi.params, snes.mu)
    candi.fitness = NaN # FIXME: use nafitness()
  end

  return snes.candidates
end

# Tell the sNES the ranking of a set of candidates.
function tell!{F}(snes::SeparableNESOpt{F}, rankedCandidates::Vector{Candidate{F}})
  u = assign_weights(rankedCandidates, snes.utilities)

  # Calc gradient
  gradient_mu = snes.last_s * u
  sq_s_minus1 = snes.last_s.^2 - 1
  gradient_sigma = sq_s_minus1 * u

  # Update the mean and sigma vectors based on the gradient
  old_mu = copy(snes.mu)
  snes.mu += snes.mu_learnrate * (snes.sigma .* gradient_mu)
  snes.sigma .*= exp(snes.sigma_learnrate / 2 * gradient_sigma)
  apply!(snes.embed, snes.mu, old_mu)

  # There is no notion of how many was better in NES so return 0
  return 0
end

# We must reorder the samples according to the order of the fitness in
# candidates!!! Or we must reorder the utilities accordingly. The latter
# is the preferred method and we can use the indices in candidates to
# accomplish it.
assign_weights{F}(candidates::Vector{Candidate{F}}, u::Vector{Float64}) =
  u[sortperm(candidates, by = fitness)]

# xNES is nice but scales badly with increasing dimension.
type XNESOpt{F,E<:EmbeddingOperator} <: NaturalEvolutionStrategyOpt
  embed::E                        # operator embedding into the search space
  lambda::Int                     # Number of samples to take per iteration
  utilities::Vector{Float64}      # Fitness utility to give to each rank
  x_learnrate::Float64
  a_learnrate::Float64
  A::Array{Float64,2}
  expA::Array{Float64,2}
  candidates::Vector{Candidate{F}}# The last sampled values, now being evaluated
  x::Individual                   # The current incumbent (aka most likely value, mu etc)
  Z::Array{Float64,2}

  function XNESOpt(embed::E; lambda::Int = 0, ini_x = nothing)
    d = numdims(search_space(embed))
    if lambda == 0
      lambda = 4 + 3*floor(Int, log(d))
    end
    if ini_x === nothing
      ini_x = rand_individual(search_space(embed))
    else
      ini_x = copy(ini_x::Individual)
    end

    new(embed, lambda, fitness_shaping_utilities_log(lambda),
      1.0, 0.5 * min(1.0 / d, 0.25),
      zeros(d, d), zeros(d, d),
      Candidate{F}[Candidate{F}(Array(Float64, d), i) for i in 1:lambda],
      ini_x, zeros(d, lambda))
  end
end

function xnes(problem::OptimizationProblem, parameters)
  params = chain(NES_DefaultOptions, parameters)
  embed = RandomBound(search_space(problem))
  XNESOpt{fitness_type(problem), typeof(embed)}(embed; lambda = params[:lambda], ini_x = params[:ini_x])
end

function ask(xnes::XNESOpt)
  copy!(xnes.expA, xnes.A)
  Base.LinAlg.expm!(xnes.expA)
  randn!(xnes.Z)
  expAZ = xnes.expA * xnes.Z

  for i in eachindex(xnes.candidates)
    candi = xnes.candidates[i]
    @inbounds for j in 1:length(candi.params)
      candi.params[j] = xnes.x[j] + expAZ[j, i]
    end
    apply!(xnes.embed, candi.params, xnes.x)
    candi.fitness = NaN # FIXME: use nafitness()
  end

  return xnes.candidates
end

function tell!{F}(xnes::XNESOpt{F}, rankedCandidates::Vector{Candidate{F}})
  u = assign_weights(rankedCandidates, xnes.utilities)

  # fixme improve memory footprint by using A_mul_B!() etc
  dA = A_mul_Bt(broadcast(*, u', xnes.Z), xnes.Z)
  for i in 1:numdims(population(xnes))
    dA[i,i] -= 1.0
  end
  dA *= xnes.a_learnrate
  xnes.A += dA

  old_x = copy(xnes.x)
  dx = xnes.x_learnrate * xnes.expA * (xnes.Z*u)
  xnes.x += dx
  apply!(xnes.embed, xnes.x, old_x)

  return 0
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
