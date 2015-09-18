abstract NaturalEvolutionStrategyOpt <: PopulationOptimizer

type SeparableNESOpt{F,E<:EmbeddingOperator} <: NaturalEvolutionStrategyOpt
  embed::E                        # operator embedding into the search space
  lambda::Int                     # Number of samples to take per iteration
  mu::Vector{Float64}             # Average position for sampling each position
  sigma::Vector{Float64}          # Average std deviation for sampling in position
  distr::Distribution             # Distribution we sample the step sizes from
  mu_learnrate::Float64
  sigma_learnrate::Float64
  max_sigma::Float64

  last_s::Array{Float64,2}        # The s values sampled in the last call to ask
  candidates::Vector{Candidate{F}}# The last sampled values, now being evaluated
  sortedUtilities::Vector{Float64}# The fitness shaping utility vector
  tmp_Utilities::Vector{Float64}   # The fitness shaping utility vector sorted by current population fitness

  function SeparableNESOpt(embed::E;
    lambda::Int = 0,
    mu_learnrate::Float64 = 1.0,
    sigma_learnrate::Float64 = 0.0,
    ini_x = nothing, max_sigma::Float64 = 1.0E+10)
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
      mu_learnrate, sigma_learnrate, max_sigma,
      zeros(d, lambda),
      Candidate{F}[Candidate{F}(Array(Float64, d), i) for i in 1:lambda],
      # Most modern NES papers use log rather than linear fitness shaping.
      fitness_shaping_utilities_log(lambda),
      @compat(Vector{Float64}(lambda)))
  end
end

population(o::NaturalEvolutionStrategyOpt) = o.candidates
numdims(o::NaturalEvolutionStrategyOpt) = numdims(search_space(o.embed))

const NES_DefaultOptions = @compat Dict{Symbol,Any}(
  :lambda => 0,              # If 0 it will be set based on the number of dimensions
  :ini_x => nothing,         # starting point, "nothing" generates random point in a search space
  :mu_learnrate => 1.0,
  :sigma_learnrate => 0.0,   # If 0.0 it will be set based on the number of dimensions
  :max_sigma => 1.0E+10       # Maximal sigma
)

function separable_nes(problem::OptimizationProblem, parameters)
  params = chain(NES_DefaultOptions, parameters)
  embed = RandomBound(search_space(problem))
  SeparableNESOpt{fitness_type(problem), typeof(embed)}(embed,
    lambda = params[:lambda],
    mu_learnrate = params[:mu_learnrate],
    sigma_learnrate = params[:sigma_learnrate],
    ini_x = params[:ini_x],
    max_sigma = params[:max_sigma])
end

calc_sigma_learnrate_for_snes(d) = (3 + log(d)) / (5 * sqrt(d))
calc_sigma_learnrate_for_nes(d) = (9 + 3 * log(d)) / (5 * d * sqrt(d))

# Get a set of new individuals to be ranked based on fitness.
function ask(snes::SeparableNESOpt)
  # Sample from N(0, 1)
  randn!(snes.last_s)

  for i in eachindex(snes.candidates)
    candi = snes.candidates[i]
    candi.index = i # reset ordering
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
  u = assign_weights!(snes.tmp_Utilities, rankedCandidates, snes.sortedUtilities)

  # Calc gradient
  gradient_mu = snes.last_s * u
  gradient_sigma = (snes.last_s.^2 - 1) * u

  # Update the mean and sigma vectors based on the gradient
  old_mu = copy(snes.mu)
  snes.mu += snes.mu_learnrate * (snes.sigma .* gradient_mu)
  snes.sigma .*= exp(snes.sigma_learnrate / 2 * gradient_sigma)
  clamp!(snes.sigma, 0.0, snes.max_sigma)
  apply!(snes.embed, snes.mu, old_mu)

  # There is no notion of how many was better in NES so return 0
  return 0
end

# given the candidates ranked by their fitness,
# the procedure returns weights sorted by the individual's index in the population
function assign_weights!{F}(weights::Vector{Float64}, rankedCandidates::Vector{Candidate{F}}, sortedWeights::Vector{Float64})
  @assert length(weights) == length(sortedWeights) && length(rankedCandidates) == length(weights)
  for i in eachindex(rankedCandidates)
    weights[rankedCandidates[i].index] = sortedWeights[i]
  end
  return weights
end

# trace current optimization state,
# Called by OptRunController trace_progress()
function trace_state(io::IO, snes::SeparableNESOpt)
    println(io, "|σ|=", norm(snes.sigma))
end

# Abstract type for a family of NES methods that represent population as
# x = μ + σ B⋅Z,
# where B is an exponential of some symmetric matrix lnB, tr(lnB)==0.0
abstract ExponentialNaturalEvolutionStrategyOpt <: NaturalEvolutionStrategyOpt

function update_candidates!(exnes::ExponentialNaturalEvolutionStrategyOpt, Z::Matrix)
  B = expm(exnes.ln_B)
  sBZ = A_mul_B!(exnes.tmp_sBZ, B, Z)
  scale!(sBZ, exnes.sigma)
  for i in eachindex(exnes.candidates)
    candi = exnes.candidates[i]
    candi.index = i # reset ordering
    @inbounds for j in 1:length(candi.params)
      candi.params[j] = exnes.x[j] + sBZ[j, i]
    end
    apply!(exnes.embed, candi.params, exnes.x)
    candi.fitness = NaN # FIXME: use nafitness()
  end

  return exnes.candidates
end

function update_parameters!(exnes::ExponentialNaturalEvolutionStrategyOpt, u::Vector)
  # TODO use syrk(Z,dA) to speed-up multiplication
  Zu = scale!(exnes.tmp_Zu, exnes.Z, u)
  ln_dB = A_mul_Bt!(exnes.tmp_lndB, Zu, exnes.Z)
  sumU = sum(u)
  dSigma = 0.0
  @inbounds for i in 1:size(ln_dB, 1)
    ln_dB[i, i] -= sumU
    dSigma += ln_dB[i, i]
  end
  dSigma /= size(ln_dB, 1)
  @inbounds for i in 1:size(ln_dB, 1)
    ln_dB[i, i] -= dSigma
  end
  scale!(ln_dB, 0.5 * exnes.B_learnrate)

  prev_x = copy!(exnes.tmp_x, exnes.x)
  dx = A_mul_B!(exnes.tmp_dx, exnes.tmp_sBZ, u)
  scale!(dx, exnes.x_learnrate)
  exnes.x += dx
  apply!(exnes.embed, exnes.x, prev_x)

  exnes.ln_B += ln_dB

  new_sigma = exnes.sigma * exp(0.5 * exnes.sigma_learnrate * dSigma)
  if isfinite(new_sigma) && new_sigma > 0
    exnes.sigma = min(new_sigma, exnes.max_sigma)
  end
  exnes
end

# identity for generic search space
ini_xnes_B(ss::SearchSpace) = eye(numdims(ss), numdims(ss))

# set B to deltas of each dimension
function ini_xnes_B(ss::RangePerDimSearchSpace)
  diag = map(log, deltas(ss))
  diag -= mean(diag)
  return full(Diagonal(diag))
end

# xNES is nice but scales badly with increasing dimension.
type XNESOpt{F,E<:EmbeddingOperator} <: ExponentialNaturalEvolutionStrategyOpt
  embed::E                        # operator embedding into the search space
  lambda::Int                     # Number of samples to take per iteration
  sortedUtilities::Vector{Float64}# Fitness utility to give to each rank
  tmp_Utilities::Vector{Float64}   # Fitness utilities assigned to current population
  x_learnrate::Float64
  sigma_learnrate::Float64
  B_learnrate::Float64
  max_sigma::Float64

  # TODO use Symmetric{Float64} to improve exponent etc calculation
  ln_B::Array{Float64,2}          # log of "covariation" matrix
  sigma::Float64                  # step size
  x::Individual                   # The current incumbent (aka most likely value, mu etc)
  Z::Array{Float64,2}             # current N(0,I) samples
  candidates::Vector{Candidate{F}}# The last sampled values, now being evaluated

  # temporary variables to minimize GC overhead
  tmp_x::Individual
  tmp_dz::Individual
  tmp_dx::Individual
  tmp_lndB::Matrix{Float64}
  tmp_Zu::Matrix{Float64}
  tmp_sBZ::Matrix{Float64}

  function XNESOpt(embed::E; lambda::Int = 0, mu_learnrate::Float64 = 1.0,
                   sigma_learnrate = 0.0, B_learnrate::Float64 = 0.0,
                   ini_x = nothing, ini_sigma::Float64 = 1.0, ini_lnB = nothing,
                   max_sigma::Float64 = 1.0E+10)
    d = numdims(search_space(embed))
    if lambda == 0
      lambda = 4 + 3*floor(Int, log(d))
    end
    if sigma_learnrate == 0.0
      sigma_learnrate = 0.5 * min(1.0 / d, 0.25) # 0.6 * (3+log(d))/(d*sqrt(d)) in other literature
    end
    if B_learnrate == 0.0
      B_learnrate = 0.6 * (3+log(d))/(d*sqrt(d)) # 0.5 * min(1.0 / d, 0.25)
    end
    if ini_x === nothing
      ini_x = rand_individual(search_space(embed))
    else
      ini_x = copy(ini_x::Individual)
      apply!(embed, ini_x, rand_individual(search_space(embed)))
    end

    new(embed, lambda, fitness_shaping_utilities_log(lambda), @compat(Vector{Float64}(lambda)),
      mu_learnrate, sigma_learnrate, B_learnrate, max_sigma,
      ini_lnB === nothing ? ini_xnes_B(search_space(embed)) : ini_lnB, ini_sigma, ini_x, zeros(d, lambda),
      Candidate{F}[Candidate{F}(Array(Float64, d), i) for i in 1:lambda],
      # temporaries
      zeros(d), zeros(d), zeros(d),
      zeros(d, d),
      zeros(d, lambda), zeros(d, lambda)
    )
  end
end

const XNES_DefaultOptions = chain(NES_DefaultOptions, @compat Dict{Symbol,Any}(
  :B_learnrate => 0.0,   # If 0.0 it will be set based on the number of dimensions
  :ini_sigma => 1.0,     # Initial sigma (step size)
  :ini_lnB => nothing    # Initial log(B)
))

function xnes(problem::OptimizationProblem, parameters)
  params = chain(XNES_DefaultOptions, parameters)
  embed = RandomBound(search_space(problem))
  XNESOpt{fitness_type(problem), typeof(embed)}(embed; lambda = params[:lambda],
                                                mu_learnrate = params[:mu_learnrate],
                                                sigma_learnrate = params[:sigma_learnrate],
                                                B_learnrate = params[:B_learnrate],
                                                ini_x = params[:ini_x],
                                                ini_sigma = params[:ini_sigma],
                                                ini_lnB = params[:ini_lnB],
                                                max_sigma = params[:max_sigma])
end

function ask(xnes::XNESOpt)
  randn!(xnes.Z)
  update_candidates!(xnes, xnes.Z)
end

function tell!{F}(xnes::XNESOpt{F}, rankedCandidates::Vector{Candidate{F}})
  u = assign_weights!(xnes.tmp_Utilities, rankedCandidates, xnes.sortedUtilities)

  update_parameters!(xnes, u)
  return 0
end

# trace current optimization state,
# Called by OptRunController trace_progress()
function trace_state(io::IO, xnes::XNESOpt)
    println(io, "sigma=", xnes.sigma,
                " |trace(ln_B)|=", trace(xnes.ln_B))
end

# Calculate the fitness shaping utilities vector using the log method.
function fitness_shaping_utilities_log(n::Int)
  u = [max(0.0, log(n / 2 + 1.0) - log(i)) for i in 1:n]
  u/sum(u) - 1/n
end

# Calculate the fitness shaping utilities vector using the steps method.
function fitness_shaping_utilities_linear(n::Int)
  # Second half has zero utility.
  treshold = floor(Int, n/2)
  second_half = zeros(n - treshold)

  # While first half's utility decreases in linear steps
  step_size = 1 / treshold
  first_half = linspace(1.0, step_size, treshold)

  # But the utilities should sum to 0, so we normalize and return
  u = vcat(first_half, second_half)
  u/sum(u) - 1/n
end
