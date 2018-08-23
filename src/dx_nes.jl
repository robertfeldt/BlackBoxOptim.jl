"""
DX-NES: distance-weighted extensions of xNES by Fukushima et al.
"""
mutable struct DXNESOpt{F,E<:EmbeddingOperator} <: ExponentialNaturalEvolutionStrategyOpt
    embed::E                        # operator embedding into the search space
    lambda::Int                     # Number of samples to take per iteration
    sortedUtilities::Vector{Float64}# Fitness utility to give to each rank
    tmp_Utilities::Vector{Float64}  # Fitness utilities assigned to current population
    x_learnrate::Float64
    B_learnrate::Float64
    sigma_learnrate::Float64
    max_sigma::Float64
    moving_threshold::Float64       # threshold of |evolutionary path| to detect movement
    evol_discount::Float64
    evol_Zscale::Float64
    uZscale::Float64                # distance scale for the weights
    evol_path::Vector{Float64}      # evolution path, discounted coordinates of population center
    # TODO use Symmetric{Float64} to inmprove exponent etc calculation
    ln_B::Matrix{Float64}           # exponential of lnB
    sigma::Float64                  # step size
    x::Individual                   # The current incumbent (aka most likely value, mu etc)
    Z::Matrix{Float64}              # current N(0,I) samples
    candidates::Vector{Candidate{F}}# The last sampled values, now being evaluated

    # temporary variables to minimize GC overhead
    tmp_x::Individual
    tmp_dz::Individual
    tmp_dx::Individual
    tmp_lndB::Matrix{Float64}
    tmp_Zu::Matrix{Float64}
    tmp_sBZ::Matrix{Float64}

    function DXNESOpt{F}(embed::E; lambda::Int = 0,
            mu_learnrate::Float64 = 1.0,
            ini_x = nothing, ini_sigma::Float64 = 1.0,
            ini_lnB = nothing,
            max_sigma::Float64 = 1.0E+10) where {F, E<:EmbeddingOperator}
        iseven(lambda) || throw(ArgumentError("lambda needs to be even"))
        d = numdims(search_space(embed))
        if lambda == 0
            lambda = 4 + 3*floor(Int, log(d))
            isodd(lambda) && (lambda += 1) # make sure lambda is even
        end
        if ini_x === nothing
            ini_x = rand_individual(search_space(embed))
        else
            ini_x = copy(ini_x::Individual)
            apply!(embed, ini_x, rand_individual(search_space(embed)))
        end
        u = fitness_shaping_utilities_log(lambda)
        moving_threshold, evol_discount, evol_Zscale = calculate_evol_path_params(d, u)

        new{F,E}(embed, lambda, u, Vector{Float64}(undef, lambda),
            mu_learnrate, 0.0, 0.0, max_sigma,
            moving_threshold, evol_discount, evol_Zscale, 0.9 + 0.15 * log(d),
            fill(NaN, d), ini_lnB === nothing ? ini_xnes_B(search_space(embed)) : ini_lnB, ini_sigma, ini_x,
            fill(NaN, d, lambda),
            [Candidate{F}(fill!(Individual(undef, d), NaN), i) for i in 1:lambda],
            # temporaries
            Vector{Float64}(undef, d), Vector{Float64}(undef, d),
            Vector{Float64}(undef, d), Matrix{Float64}(undef, d, d),
            Matrix{Float64}(undef, d, lambda), Matrix{Float64}(undef, d, lambda)
        )
    end
end

function trace_state(io::IO, dxnes::DXNESOpt, mode::Symbol)
    evol_path_norm = norm(dxnes.evol_path)
    println(io,
            "σ=", dxnes.sigma,
            " η[x]=", dxnes.x_learnrate,
            " η[σ]=", dxnes.sigma_learnrate,
            " η[B]=", dxnes.B_learnrate,
            " |tr(ln_B)|=", abs(tr(dxnes.ln_B)),
            " |path|=", evol_path_norm,
            " speed=", evol_path_norm/dxnes.moving_threshold)
end

const DXNES_DefaultOptions = chain(NES_DefaultOptions, ParamsDict(
    :ini_sigma => 1.0,      # Initial sigma (step size)
    :ini_lnB => nothing     # Initial log(B) (log of parameters covariation)
))

function dxnes(problem::OptimizationProblem, parameters)
    params = chain(DXNES_DefaultOptions, parameters)
    embed = RandomBound(search_space(problem))
    DXNESOpt{fitness_type(problem)}(embed; lambda = params[:lambda],
                                    ini_x = params[:ini_x],
                                    ini_sigma = params[:ini_sigma],
                                    ini_lnB = params[:ini_lnB],
                                    max_sigma = params[:max_sigma])
end

function ask(dxnes::DXNESOpt)
    hlambda = dxnes.lambda÷2
    @inbounds for i in 1:hlambda
        for j in 1:size(dxnes.Z, 1)
            dxnes.Z[j, i] = randn()
            dxnes.Z[j, i+hlambda] = -dxnes.Z[j, i]
        end
    end
    update_candidates!(dxnes, dxnes.Z)
end

function tell!(dxnes::DXNESOpt{F}, rankedCandidates::Vector{Candidate{F}}) where F
    u = assign_weights!(dxnes.tmp_Utilities, rankedCandidates, dxnes.sortedUtilities)
    dxnes.evol_path *= dxnes.evol_discount
    # We'll take the small perf hit for now just so this can run also on pre rc2 julia 0.7s
    dxnes.evol_path += dxnes.evol_Zscale * dropdims(wsum(dxnes.Z, u, 2), dims=2)
    evol_speed = norm(dxnes.evol_path)/dxnes.moving_threshold
    if evol_speed > 1.0
        # the center is moving, adjust weights
        u = assign_distance_weights!(dxnes.tmp_Utilities, dxnes.uZscale, rankedCandidates, dxnes.Z)
    end
    update_learning_rates!(dxnes, evol_speed)
    update_parameters!(dxnes, u)

    return 0
end

"""
Calculate the parameters for evolutionary path.

Returns the tuple:
  * moving threshold
  * path discount
  * current Z scale
"""
function calculate_evol_path_params(n::Int, u::Vector{Float64})
    lambda = length(u)
    mu = 1/sum(i -> (u[i]+1/lambda)^2, 1:(lambda÷2))
    c = (mu + 2.0)/(sqrt(n)*(n+mu+5.0))
    discount = 1-c
    Zscale = sqrt(c*(2-c)*mu)
    # Expected length of evol.path, if uZ ~ α N(0,I), based on AR model:
    # E|path|^2 = α*Zscale^2 / (1-discount^2),
    # where α ≈ 1 if the population does not cover the minimum and lies on the slope
    # (so that |uZ| is large), but σ is too large and the method is diverging
    threshold = Zscale*sqrt(n/(1-discount^2))
    #info("DX-NES params: moving_threshold=", threshold,
    #     " μ=", mu, " c=", c, " discount=", discount, " Zscale=", Zscale)
    return threshold, discount, Zscale
end

function assign_distance_weights!(weights::Vector{Float64}, scale::Float64,
                                  rankedCandidates::Vector{Candidate{F}},
                                  Z::Matrix{Float64}) where F
    @assert length(weights) == length(rankedCandidates) == size(Z, 2)
    lambda = size(Z, 2)
    u0 = log(0.5*lambda+1)
    for i in 1:lambda
        cand_ix = rankedCandidates[i].index
        weights[cand_ix] = max(u0-log(i), 0) * exp(scale*norm(view(Z, cand_ix)))
    end
    weights .= weights./sum(weights) .- 1/lambda
    return weights
end

function update_learning_rates!(dxnes::DXNESOpt, evol_speed::Float64)
    n = numdims(dxnes)
    if evol_speed > 1.0
        dxnes.sigma_learnrate = 1.0
        dxnes.B_learnrate = (dxnes.lambda + 2n)/(dxnes.lambda+2*n^2+100) *
                            min(1.0, sqrt(dxnes.lambda/numdims(dxnes)))
    else
        dxnes.sigma_learnrate = 1.0 + dxnes.lambda/(dxnes.lambda+2*n)
        if evol_speed > 0.1
            dxnes.sigma_learnrate *= 0.5
        end
        dxnes.B_learnrate = dxnes.lambda/(dxnes.lambda+2*n^2+100)
    end
    return dxnes
end
