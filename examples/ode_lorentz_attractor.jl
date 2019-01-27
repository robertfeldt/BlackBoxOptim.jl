# Let's try parameter estimation for an ODE, here the famous Lorentz attractor.
# We base our implementation on the code from Paulo Marques available at:
#   https://github.com/pjpmarques/Julia-Modeling-the-World/blob/master/Lorenz%20Attractor.ipynb

function lorentz_ode!(dr::AbstractVector{Float64},
                      params::AbstractVector{Float64},
                      r::AbstractVector{Float64})
    @assert length(dr) == length(r) == 3

    # Get individual params for this ODE
    @inbounds sigma, rho, beta = params

    # Extract the coordinates from the r vector
    @inbounds x, y, z = r

    # The Lorenz equations, put the derivative into dr
    @inbounds dr[1] = sigma*(y - x)       # dx/dt
    @inbounds dr[2] = x*(rho - z) - y     # dy/dt
    @inbounds dr[3] = x*y - beta*z        # dz/dt

    return dr
end

# Define time vector and interval grid.
t = range(0.0, step=0.001, stop=100.0)
# But in order to compare to [Xiang2015] paper we use their choice instead:
t_Xiang2015 = range(0.0, step=0.01, length=301)

# Initial position in space
r0 = [0.1, 0.0, 0.0]

# Constants sigma, rho and beta. In a real ODE problem these would not be known
# and would be estimated.
sigma = 10.0
rho   = 28.0
beta  = 8.0/3.0
real_params = [sigma, rho, beta]

# "Play" an ODE from a starting point and into the future given a sequence of time steps.
# Put the restults into `states`
function calc_states!(states::AbstractMatrix{Float64},
    params::AbstractVector{Float64}, odefunc!::Function,
    startx::AbstractVector{Float64}, times::AbstractVector{Float64})

    numsamples = length(times)
    @assert size(states) == (length(startx), numsamples)

    states[:, 1] = startx
    tprev = times[1]
    deriv = similar(startx)
    for i in 2:numsamples
        tnow = times[i]
        stateprev = view(states, :, i-1)
        odefunc!(deriv, params, stateprev)
        @inbounds states[:, i] .= stateprev .+ deriv .* (tnow - tprev)
        tprev = tnow
    end

    return states
end

calc_states(params::AbstractVector{Float64}, odefunc!::Function,
    startx::AbstractVector{Float64}, times::AbstractVector{Float64}) =
    calc_states!(Matrix{Float64}(undef, length(startx), length(times)),
                 params, odefunc!, startx, times)

# RSS = Residual Sum of Squares, columnwise
rss(actual::AbstractMatrix{Float64}, estimated::AbstractMatrix{Float64}) =
    @inbounds sum(i -> abs2(actual[i] - estimated[i]), eachindex(actual, estimated))

# Calculate the actual/original state vectors that we will use for parameter
# estimation:
origstates = calc_states(real_params, lorentz_ode!, r0, t)
origstates_Xiang2015 = calc_states(real_params, lorentz_ode!, r0, t_Xiang2015)

function subsample(origstates::AbstractMatrix{Float64},
                   times::AbstractVector{Float64};
                   lenratio = 0.25)
    N = length(times)
    @assert size(origstates, 2) == N
    stopidx = round(Int, lenratio*N)
    indexes = 1:stopidx
    return origstates[:, indexes], times[indexes]
end

# The [Xiang2015] paper, https://www.hindawi.com/journals/ddns/2015/740721/,
# used these param bounds:
Xiang2015Bounds = [(9., 11.), (20., 30.), (2., 3.)]

# Now we can optimize using BlackBoxOptim
using BlackBoxOptim

# store temporary states for fitness calculation
const tmpstates_pool = Dict{NTuple{2, Int}, Vector{Matrix{Float64}}}()

# get the RSS between the origstates trajectory and ODE solution for given params
function lorentz_fitness(params::AbstractVector{Float64},
                         origstates::AbstractMatrix{Float64},
                         times::AbstractVector{Float64})
    # get the state matrix from the pool of the proper size
    statesize = size(origstates)
    states_pool = get!(() -> Vector{Matrix{Float64}}(), tmpstates_pool, statesize)
    states = !isempty(states_pool) ? pop!(states_pool) : Matrix{Float64}(undef, statesize...)
    # solve ODE for given params
    calc_states!(states, params, lorentz_ode!, r0, times)
    push!(states_pool, states) # return states matrix back to the pool
    return rss(origstates, states)
end

# But optimizing all states in each optimization step is too much so lets
# sample a small subset and use for first opt iteration.
origstates1, times1 = subsample(origstates, t; lenratio=0.04); # Sample only first 4%
origstates1_Xiang2015, times1_Xiang2015 = subsample(origstates_Xiang2015, t_Xiang2015; lenratio=1.00);

res1 = bboptimize(params -> lorentz_fitness(params, origstates1, times1);
    SearchRange = Xiang2015Bounds, MaxSteps = 8e3)

res2 = bboptimize(params -> lorentz_fitness(params, origstates1_Xiang2015, times1_Xiang2015);
    SearchRange = Xiang2015Bounds, MaxSteps = 11e3) # They allow 12k fitness evals for 3-param estimation

# But lets also try with less tight bounds
LooserBounds = [(0., 22.), (0., 60.), (1., 6.)]
res3 = bboptimize(params -> lorentz_fitness(params, origstates1_Xiang2015, times1_Xiang2015);
    SearchRange = LooserBounds, MaxSteps = 11e3) # They allow 12k fitness evals for 3-param estimation

@info "Results on the long time sequence from Paulo Marques:"
estfitness = lorentz_fitness(best_candidate(res1), origstates, t)
@show estfitness best_candidate(res1) best_fitness(res1)
origfitness = lorentz_fitness(real_params, origstates, t)
@show origfitness real_params

@info "Results on the short time sequence used in [Xiang2015] paper:"
estfitness = lorentz_fitness(best_candidate(res2), origstates_Xiang2015, t_Xiang2015)
@show estfitness best_candidate(res2) best_fitness(res2)
origfitness = lorentz_fitness(real_params, origstates_Xiang2015, t_Xiang2015)
@show origfitness real_params

@info "Results on the short time sequence used in [Xiang2015] paper, but with looser bounds:"
estfitness = lorentz_fitness(best_candidate(res3), origstates_Xiang2015, t_Xiang2015)
@show estfitness best_candidate(res3) best_fitness(res3)
origfitness = lorentz_fitness(real_params, origstates_Xiang2015, t_Xiang2015)
@show origfitness real_params
