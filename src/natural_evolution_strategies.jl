abstract type NaturalEvolutionStrategyOpt <: PopulationOptimizer end

"""
Separable Natural Evolution Strategy (sNES) optimizer.
"""
mutable struct SeparableNESOpt{F,E<:EmbeddingOperator} <: NaturalEvolutionStrategyOpt
    embed::E                        # operator embedding into the search space
    lambda::Int                     # number of samples to take per iteration
    mu::Vector{Float64}             # center of the population
    sigma::Vector{Float64}          # average std deviation for sampling in position
    distr::Distribution             # distribution to sample the step sizes from
    mu_learnrate::Float64
    sigma_learnrate::Float64
    max_sigma::Float64

    last_s::Matrix{Float64}         # `s` values sampled in the last call to `ask()`
    candidates::Vector{Candidate{F}}# the last sampled values, now being evaluated
    sortedUtilities::Vector{Float64}# the fitness shaping utility vector
    tmp_Utilities::Vector{Float64}   # the fitness shaping utility vector sorted by current population fitness

    function SeparableNESOpt{F}(
            embed::E;
            lambda::Int = 0,
            mu_learnrate::Float64 = 1.0,
            sigma_learnrate::Float64 = 0.0,
            ini_x = nothing, max_sigma::Float64 = 1.0E+10
    ) where {F, E<:EmbeddingOperator}
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

        new{F,E}(embed, lambda,
                 ini_x, fill(1.0, d),
                 Normal(0, 1),
                 mu_learnrate, sigma_learnrate, max_sigma,
                 fill(0.0, d, lambda),
                 [Candidate{F}(fill!(Individual(undef, d), NaN), i) for i in 1:lambda],
                 # Most modern NES papers use log rather than linear fitness shaping.
                 fitness_shaping_utilities_log(lambda),
                 Vector{Float64}(undef, lambda))
    end
end

population(o::NaturalEvolutionStrategyOpt) = o.candidates
numdims(o::NaturalEvolutionStrategyOpt) = numdims(search_space(o.embed))

const NES_DefaultOptions = ParamsDict(
    :lambda => 0,              # If 0 it will be set based on the number of dimensions
    :ini_x => nothing,         # starting point, "nothing" generates random point in a search space
    :mu_learnrate => 1.0,
    :sigma_learnrate => 0.0,   # If 0.0 it will be set based on the number of dimensions
    :max_sigma => 1.0E+10      # Maximal sigma
)

function separable_nes(problem::OptimizationProblem, parameters)
    params = chain(NES_DefaultOptions, parameters)
    embed = RandomBound(search_space(problem))
    SeparableNESOpt{fitness_type(problem)}(embed,
        lambda = params[:lambda],
        mu_learnrate = params[:mu_learnrate],
        sigma_learnrate = params[:sigma_learnrate],
        ini_x = params[:ini_x],
        max_sigma = params[:max_sigma])
end

calc_sigma_learnrate_for_snes(d) = (3 + log(d)) / (5 * sqrt(d))
calc_sigma_learnrate_for_nes(d) = (9 + 3 * log(d)) / (5 * d * sqrt(d))

function ask(snes::SeparableNESOpt)
    # Sample from N(0, 1)
    randn!(snes.last_s)

    for i in eachindex(snes.candidates)
        candi = snes.candidates[i]
        candi.index = i # reset ordering
        @inbounds for j in eachindex(candi.params)
            candi.params[j] = snes.mu[j] + snes.sigma[j] * snes.last_s[j, i]
        end
        apply!(snes.embed, candi.params, snes.mu)
        candi.fitness = NaN # FIXME: use nafitness()
    end

    return snes.candidates
end

function tell!(snes::SeparableNESOpt{F}, rankedCandidates::AbstractVector{<:Candidate{F}}) where F
    u = assign_weights!(snes.tmp_Utilities, rankedCandidates, snes.sortedUtilities)

    # Calc gradient
    gradient_mu = snes.last_s * u
    gradient_sigma = (abs2.(snes.last_s) .- 1) * u

    # Update the mean and sigma vectors based on the gradient
    old_mu = copy(snes.mu)
    snes.mu += snes.mu_learnrate * (snes.sigma .* gradient_mu)
    @inbounds for i in eachindex(snes.sigma)
        snes.sigma[i] = clamp(snes.sigma[i] * exp(0.5 * snes.sigma_learnrate * gradient_sigma[i]),
                              0.0, snes.max_sigma)
    end
    apply!(snes.embed, snes.mu, old_mu)

    # There is no notion of how many was better in NES so return 0
    return 0
end

"""
    assign_weights!(weights, rankedCandidates, sortedWeights)

Assigns the candidate `weights` according to the candidate index.
`rankedCandidates` are ranked by their fitness, `sortedWeights` are
the corresponding weights.

Returns candidate weights sorted by the individual's index in the population.
"""
function assign_weights!(weights::Vector{Float64},
                         rankedCandidates::AbstractVector{<:Candidate},
                         sortedWeights::Vector{Float64})
    @assert length(weights) == length(sortedWeights) && length(rankedCandidates) == length(weights)
    for i in eachindex(rankedCandidates)
        weights[rankedCandidates[i].index] = sortedWeights[i]
    end
    return weights
end

#  Traces the current `sNES` optimization state.
#  Called by `OptRunController` `trace_progress()`.
function trace_state(io::IO, snes::SeparableNESOpt, mode::Symbol)
    println(io, "|σ|=", norm(snes.sigma))
end

"""
Abstract type for a family of NES methods that represent population as
```
x = μ + σ B⋅Z,
```
where `B` is an exponential of some symmetric matrix `lnB`, `tr(lnB)==0.0`
"""
abstract type ExponentialNaturalEvolutionStrategyOpt <: NaturalEvolutionStrategyOpt end

function update_candidates!(exnes::ExponentialNaturalEvolutionStrategyOpt, Z::Matrix)
    B = exp(exnes.ln_B)
    sBZ = mul!(exnes.tmp_sBZ, B, Z)
    sBZ .*= exnes.sigma
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

function update_parameters!(exnes::ExponentialNaturalEvolutionStrategyOpt, u::AbstractVector)
    # TODO use syrk(Z,dA) to speed-up multiplication
    Zu = broadcast!(*, exnes.tmp_Zu, exnes.Z, transpose(u))
    ln_dB = mul!(exnes.tmp_lndB, Zu, transpose(exnes.Z))
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
    ln_dB .*= 0.5 * exnes.B_learnrate

    prev_x = copyto!(exnes.tmp_x, exnes.x)
    dx = mul!(exnes.tmp_dx, exnes.tmp_sBZ, u)
    dx .*= exnes.x_learnrate
    exnes.x += dx
    apply!(exnes.embed, exnes.x, prev_x)

    exnes.ln_B += ln_dB

    new_sigma = exnes.sigma * exp(0.5 * exnes.sigma_learnrate * dSigma)
    if isfinite(new_sigma) && new_sigma > 0
        exnes.sigma = min(new_sigma, exnes.max_sigma)
    end
    return exnes
end

" Identity for generic search space "
ini_xnes_B(ss::SearchSpace) = Matrix{Float64}(I, numdims(ss), numdims(ss))

"""
Calculates the initial ``log B`` matrix for `xNES` based on the deltas of each dimension.
"""
function ini_xnes_B(ss::RangePerDimSearchSpace)
    diag = log.(deltas(ss))
    diag .-= mean(diag)
    return Matrix{Float64}(Diagonal(diag))
end

"""
`xNES` method.

Nice but scales badly with increasing dimensions.
"""
mutable struct XNESOpt{F,E<:EmbeddingOperator} <: ExponentialNaturalEvolutionStrategyOpt
    embed::E                        # operator embedding into the search space
    lambda::Int                     # number of samples to take per iteration
    sortedUtilities::Vector{Float64}# the fitness shaping utility vector
    tmp_Utilities::Vector{Float64}   # the fitness shaping utility vector sorted by current population fitness
    x_learnrate::Float64
    sigma_learnrate::Float64
    B_learnrate::Float64
    max_sigma::Float64

    # TODO use Symmetric{Float64} to improve exponent etc calculation
    ln_B::Matrix{Float64}           # `log` of "covariation" matrix
    sigma::Float64                  # step size
    x::Individual                   # center of the population (aka most likely value, `mu` etc)
    Z::Matrix{Float64}              # current `N(0,I)` samples
    candidates::Vector{Candidate{F}}# the last sampled values, now being evaluated

    # temporary variables to minimize GC overhead
    tmp_x::Individual
    tmp_dz::Individual
    tmp_dx::Individual
    tmp_lndB::Matrix{Float64}
    tmp_Zu::Matrix{Float64}
    tmp_sBZ::Matrix{Float64}

    function XNESOpt{F}(embed::E; lambda::Int = 0, mu_learnrate::Float64 = 1.0,
            sigma_learnrate = 0.0, B_learnrate::Float64 = 0.0,
            ini_x = nothing, ini_sigma::Float64 = 1.0, ini_lnB = nothing,
            max_sigma::Float64 = 1.0E+10) where {F, E<:EmbeddingOperator}
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

        new{F,E}(embed, lambda, fitness_shaping_utilities_log(lambda),
                 fill!(Vector{Float64}(undef, lambda), NaN),
                 mu_learnrate, sigma_learnrate, B_learnrate, max_sigma,
                 ini_lnB === nothing ? ini_xnes_B(search_space(embed)) : ini_lnB, ini_sigma, ini_x, fill(0.0, d, lambda),
                 [Candidate{F}(fill!(Individual(undef, d), NaN), i) for i in 1:lambda],
                  # temporaries
                  Vector{Float64}(undef, d), Vector{Float64}(undef, d), Vector{Float64}(undef, d),
                  Matrix{Float64}(undef, d, d),
                  Matrix{Float64}(undef, d, lambda), Matrix{Float64}(undef, d, lambda),
        )
    end
end

const XNES_DefaultOptions = chain(NES_DefaultOptions, ParamsDict(
    :B_learnrate => 0.0,   # If 0.0 it will be set based on the number of dimensions
    :ini_sigma => 1.0,     # Initial sigma (step size)
    :ini_lnB => nothing    # Initial log(B)
))

function xnes(problem::OptimizationProblem, parameters)
    params = chain(XNES_DefaultOptions, parameters)
    embed = RandomBound(search_space(problem))
    F = fitness_type(problem)
    XNESOpt{F}(embed; lambda = params[:lambda],
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

function tell!(xnes::XNESOpt{F}, rankedCandidates::Vector{Candidate{F}}) where F
    u = assign_weights!(xnes.tmp_Utilities, rankedCandidates, xnes.sortedUtilities)

    update_parameters!(xnes, u)
    return 0
end

function trace_state(io::IO, xnes::XNESOpt, mode::Symbol)
    println(io, "sigma=", xnes.sigma,
                " |trace(ln_B)|=", tr(xnes.ln_B))
end

"""
    fitness_shaping_utilities_log(n)

Calculate the `n`-dimensional fitness shaping utilities vector using the "log" method.
"""
function fitness_shaping_utilities_log(n::Int)
    u = max.(0.0, log(n / 2 + 1.0) .- log.(1:n))
    return u./sum(u) .- 1/n
end

"""
    fitness_shaping_utilities_linear(n)

Calculate the `n`-dimensional fitness shaping utilities vector
using the "steps" method.
"""
function fitness_shaping_utilities_linear(n::Int)
    # Second half has zero utility.
    threshold = n÷2
    second_half = fill(0.0, n - threshold)

    # While first half's utility decreases in linear steps
    first_half = range(1.0, stop=1/threshold, length=threshold)

    # But the utilities should sum to 0, so we normalize and return
    u = vcat(first_half, second_half)
    return u./sum(u) .- 1/n
end
