# Let's try parameter estimation for an ODE, here the famous Lorentz attractor.
# We base our implementation on the code from Paulo Marques available at:
#   https://github.com/pjpmarques/Julia-Modeling-the-World/blob/master/Lorenz%20Attractor.ipynb
function lorentz_equations(params::Vector{Float64}, r::Vector{Float64})
    # Get individual params for this ODE
    sigma, rho, beta = params

    # Extract the coordinates from the r vector
    x, y, z = r
    
    # The Lorenz equations
    dx_dt = sigma*(y - x)
    dy_dt = x*(rho - z) - y
    dz_dt = x*y - beta*z
    
    # Return the derivatives as a vector
    Float64[dx_dt, dy_dt, dz_dt]
end

# Define time vector and interval grid.
dt = 0.001
tf = 100.0
tinterval = 0:dt:tf
t  = collect(tinterval)

# But in order to compare to [Xiang2015] paper we use their choice instead:
h = 0.01
M = 300
tstart = 0.0
tstop = tstart + M * h
tinterval_Xiang2015 = 0:h:tstop
t_Xiang2015 = collect(tinterval_Xiang2015)

# Initial position in space
r0 = [0.1; 0.0; 0.0]

# Constants sigma, rho and beta. In a real ODE problem these would not be known 
# and would be estimated.
sigma = 10.0
rho   = 28.0
beta  = 8.0/3.0
real_params = [sigma, rho, beta]

# "Play" an ODE from a starting point and into the future given a sequence of time steps.
function calc_state_vectors(params::Vector{Float64}, odefunc::Function, 
    startx::Vector{Float64}, times::Vector{Float64}; states = nothing)

    numsamples = length(times)
    if states == nothing
        states = Array(Float64, length(startx), numsamples)
    end

    states[:, 1] = startx
    tprev = times[1]
    for i in 2:numsamples
        @inbounds derivatives = odefunc(params, states[:, (i-1)])
        @inbounds tnow = times[i]
        @inbounds states[:, i] = states[:, (i-1)] .+ derivatives * (tnow - tprev)
        tprev = tnow
    end

    return states
end

# RSS = Residual Sum of Squares, columnwise
function rss(actual::Array{Float64, 2}, estimated::Array{Float64, 2})
    M = size(actual, 2)
    sumsq = 0.0
    for i in 1:M
        @inbounds sumsq += sumabs2(actual[:, i] .- estimated[:, i])
    end
    sumsq
end

# Calculate the actual/original state vectors that we will use for parameter
# estimation:
origstates = calc_state_vectors(real_params, lorentz_equations, r0, t)
origstates_Xiang2015 = calc_state_vectors(real_params, lorentz_equations, r0, t_Xiang2015)

function subsample(origstates::Array{Float64, 2}, times::Vector{Float64}, lenratio = 0.25)
    @assert size(origstates, 2) == length(times)
    N = length(times)
    stopidx = round(Int, lenratio*N)
    indexes = 1:stopidx
    return origstates[:, indexes], times[indexes]
end

# The [Xiang2015] paper, https://www.hindawi.com/journals/ddns/2015/740721/, 
# used these param bounds:
Xiang2015Bounds = Tuple{Float64, Float64}[(9, 11), (20, 30), (2, 3)]

# Now we can optimize using BlackBoxOptim
using BlackBoxOptim

function lorentz_fitness(params::Vector{Float64}, origstates::Array{Float64, 2}, times::Vector{Float64})
    states = calc_state_vectors(params, lorentz_equations, r0, times)
    return rss(origstates, states)
end

# But optimizing all states in each optimization step is too much so lets
# sample a small subset and use for first opt iteration.
origstates1, times1 = subsample(origstates, t, 0.04); # Sample only first 4%
origstates1_Xiang2015, times1_Xiang2015 = subsample(origstates_Xiang2015, t_Xiang2015, 1.00);

res1 = bboptimize(params -> lorentz_fitness(params, origstates1, times1); 
    SearchRange = Xiang2015Bounds, MaxSteps = 8e3)

res2 = bboptimize(params -> lorentz_fitness(params, origstates1_Xiang2015, times1_Xiang2015); 
    SearchRange = Xiang2015Bounds, MaxSteps = 11e3) # They allow 12k fitness evals for 3-param estimation

# But lets also try with less tight bounds
LooserBounds = Tuple{Float64, Float64}[(0, 22), (0, 60), (1, 6)]
res3 = bboptimize(params -> lorentz_fitness(params, origstates1_Xiang2015, times1_Xiang2015); 
    SearchRange = LooserBounds, MaxSteps = 11e3) # They allow 12k fitness evals for 3-param estimation

println("Results on the long time sequence from Paulo Marques:")
estfitness = lorentz_fitness(best_candidate(res1), origstates, t)
@show (estfitness, best_candidate(res1), best_fitness(res1))
origfitness = lorentz_fitness(real_params, origstates, t)
@show (origfitness, real_params)

println("Results on the short time sequence used in [Xiang2015] paper:")
estfitness = lorentz_fitness(best_candidate(res2), origstates_Xiang2015, t_Xiang2015)
@show (estfitness, best_candidate(res2), best_fitness(res2))
origfitness = lorentz_fitness(real_params, origstates_Xiang2015, t_Xiang2015)
@show (origfitness, real_params)

println("Results on the short time sequence used in [Xiang2015] paper, but with looser bounds:")
estfitness = lorentz_fitness(best_candidate(res3), origstates_Xiang2015, t_Xiang2015)
@show (estfitness, best_candidate(res3), best_fitness(res3))
origfitness = lorentz_fitness(real_params, origstates_Xiang2015, t_Xiang2015)
@show (origfitness, real_params)