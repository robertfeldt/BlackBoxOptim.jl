# DX-NES: distance-weighted extensions of xNES by Fukushima et al
type DXNESOpt{F,E<:EmbeddingOperator} <: ExponentialNaturalEvolutionStrategyOpt
  embed::E                        # operator embedding into the search space
  lambda::Int                     # Number of samples to take per iteration
  sortedUtilities::Vector{Float64}# Fitness utility to give to each rank
  tmp_Utilities::Vector{Float64}   # Fitness utilities assigned to current population
  x_learnrate::Float64
  B_learnrate::Float64
  sigma_learnrate::Float64
  moving_threshold::Float64       # threshold of evolutionary path movement
  evol_discount::Float64
  evol_Zscale::Float64
  evol_path::Vector{Float64}      # evolution path, discounted coordinates of population center
  # TODO use Symmetric{Float64} to inmprove exponent etc calculation
  ln_B::Array{Float64,2}          # exponential of lnB
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

  function DXNESOpt(embed::E; lambda::Int = 0,
                   mu_learnrate::Float64 = 1.0,
                   ini_x = nothing, ini_sigma::Float64 = 1.0)
    iseven(lambda) || throw(ArgumentError("lambda needs to be even"))
    d = numdims(search_space(embed))
    if lambda == 0
      lambda = 4 + 3*floor(Int, log(d))
    end
    if ini_x === nothing
      ini_x = rand_individual(search_space(embed))
    else
      ini_x = copy(ini_x::Individual)
    end
    u = fitness_shaping_utilities_log(lambda)
    evol_discount, evol_Zscale = calculate_evol_path_params(d, u)

    new(embed, lambda, u, @compat(Vector{Float64}(lambda)),
      mu_learnrate, 0.0, 0.0,
      mean(Chi(d)), evol_discount, evol_Zscale, zeros(d),
      ini_xnes_B(search_space(embed)), ini_sigma, ini_x,
      zeros(d, lambda),
      Candidate{F}[Candidate{F}(Array(Float64, d), i) for i in 1:lambda],
      # temporaries
      zeros(d), zeros(d), zeros(d),
      zeros(d, d),
      zeros(d, lambda), zeros(d, lambda)
    )
  end
end

# trace current optimization state,
# Called by OptRunController trace_progress()
function trace_state(io::IO, dxnes::DXNESOpt)
    println(io,
            "σ=", dxnes.sigma,
            " η[x]=", dxnes.x_learnrate,
            " η[σ]=", dxnes.sigma_learnrate,
            " η[B]=", dxnes.B_learnrate,
            " |tr(ln_B)|=", abs(trace(dxnes.ln_B)),
            " |path|=", norm(dxnes.evol_path),
            " (", is_moving(dxnes) ? "moving" : "stopped", ")")
end

const DXNES_DefaultOptions = chain(NES_DefaultOptions, @compat Dict{Symbol,Any}(
  :ini_sigma => 1.0      # Initial sigma (step size)
))

function dxnes(problem::OptimizationProblem, parameters)
  params = chain(DXNES_DefaultOptions, parameters)
  embed = RandomBound(search_space(problem))
  DXNESOpt{fitness_type(problem), typeof(embed)}(embed; lambda = params[:lambda],
                                                 ini_x = params[:ini_x],
                                                 ini_sigma = params[:ini_sigma])
end

function ask(dxnes::DXNESOpt)
  hlambda = dxnes.lambda÷2
  for i in 1:hlambda
    for j in 1:size(dxnes.Z, 1)
      dxnes.Z[j, i] = randn()
      dxnes.Z[j, i+hlambda] = -dxnes.Z[j, i]
    end
  end
  update_candidates!(dxnes, dxnes.Z)
end

function tell!{F}(dxnes::DXNESOpt{F}, rankedCandidates::Vector{Candidate{F}})
  u = assign_weights!(dxnes.tmp_Utilities, rankedCandidates, dxnes.sortedUtilities)
  dxnes.evol_path *= dxnes.evol_discount
  dxnes.evol_path += dxnes.evol_Zscale * squeeze(wsum(dxnes.Z, u, 2), 2)
  if is_moving(dxnes)
    # the center is moving, adjust weights
    u = assign_distance_weights!(dxnes.tmp_Utilities, rankedCandidates, dxnes.Z)
  end
  update_learning_rates!(dxnes)
  update_parameters!(dxnes, u)

  return 0
end

is_moving(dxnes::DXNESOpt) = norm(dxnes.evol_path) > dxnes.moving_threshold

function calculate_evol_path_params(n::Int64, u::Vector{Float64})
  lambda = length(u)
  mu = 1/sum(i -> (u[i]+1/lambda)^2, 1:(lambda÷2))
  c = (mu + 2.0)/(sqrt(n)*(n+mu+5.0))
  return (1-c), sqrt(c*(2-c)*mu)
end

function assign_distance_weights!{F}(weights::Vector{Float64}, rankedCandidates::Vector{Candidate{F}}, Z::Matrix{Float64})
  @assert length(weights) == length(rankedCandidates) == size(Z, 2)
  lambda = size(Z, 2)
  u0 = log(0.5*lambda+1)
  for i in 1:lambda
    cand_ix = rankedCandidates[i].index
    weights[cand_ix] = max(u0-log(i), 0) * exp(norm(Z[:,cand_ix]))
  end
  wsum = sum(weights)
  for i in 1:lambda
    weights[i] = weights[i]/wsum - 1/lambda
  end
  return weights
end

function update_learning_rates!(dxnes::DXNESOpt)
  n = numdims(dxnes)
  if is_moving(dxnes)
    dxnes.sigma_learnrate = 1.0
    dxnes.B_learnrate = (dxnes.lambda + 2n)/(dxnes.lambda+2*n^2+100) *
                        min(1, sqrt(dxnes.lambda/numdims(dxnes)))
  else
    dxnes.sigma_learnrate = 1.0 + dxnes.lambda/(dxnes.lambda+2*n)
    if norm(dxnes.evol_path) > 0.1*dxnes.moving_threshold
      dxnes.sigma_learnrate *= 0.5
    end
    dxnes.B_learnrate = dxnes.lambda/(dxnes.lambda+2*n^2+100)
  end
  dxnes
end
