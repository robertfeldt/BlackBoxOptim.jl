# This is an implementation of the CMSA-ES as described in the paper:
#  H. Beyer and B. Sendhoff, "Covariance Matrix Adaptation Revisited – 
#  the CMSA Evolution Strategy", 2008
#

#type CMSA_ES <: PopulationOptimizer
#  decomposer::Function    # Function that decomposes covar matrix into sqrtC
#
#  # Constants that we save for performance reasons
#  n::Int64
#  mu::Int64
#  lambda::Int64
#  tau::Float64
#  tau_c::Float64
#
#  # State
#  sigma::Float64            # Current mutation strength, sigma
#  sigmas::Array{Float64,2}  # Sigma per individual in last generated population
#  C::Array{Float64,2}       # Covariance matrix
#  sqrtC::Array{Float64,2}   # sqrtC' * sqrtC == C, i.e. some composition
#
#  CMSA_ES(n, mu, lambda, decomposeFunc = chol) = begin
#
#    tau = 1 / sqrt(2*n)           # Equation (1) on page 5 in Beyer2008 
#    tau_c = 1 + n * (n + 1) / 2   # Equation (2) on page 5 in Beyer2008
#
#    sigma = 1.0                   # Mutation strength, sigma (self-adapted) 
#    sigmas = sigma * ones(1, n)
#    C = eye(n, n)                 # Covariance matrix
#
#    new(n, mu, lambda, decomposeFunc, tau, tau_c, 
#      sigma, sigmas, C, decomposeFunc(C))
#
#  end
#end
#
#
#
## Get a set of new individuals to be ranked based on fitness.
#function ask(o::CMSA_ES)
#  o.sigmas = o.sigma * exp(o.tau * randn(1, o.lambda))
#  o.s = randn(o.n, o.n) * o.sqrtC
#end

function normalize_utilities(utilities)
  utilities / sum(utilities)
end

function linear_utilities(mu, lambda)
  normalize_utilities(hcat(ones(1, mu), zeros(1, lambda-mu)))
end

function log_utilities(mu, lambda)
  normalize_utilities(hcat((log(mu+1) - log(1:mu))', zeros(1, lambda-mu)))
end

# Optimize the n-dimensional objective func with a (mu,lambda) CMSA-ES.
function cmsa_es(n, func; 
  mu = false, lambda = false, num_reps = 5e2,
  covarMatrixSampler = EigenCovarSampler,
  utilitiesFunc = linear_utilities)

  if lambda == false
    # In the paper they try population sizes 8, 4n and 4*n*n, lets try
    # in between:
    lambda = int(8 + (4*n*n - 8) * rand())
    # We should probably use an IPOP scheme for restarts...
  end

  if mu == false
    mu = max(int(sqrt(lambda)), 2)
  end

  tau = 1 / sqrt(2*n)           # Equation (1) on page 5 in Beyer2008 
  tau_c = 1 + n * (n + 1) / 2   # Equation (2) on page 5 in Beyer2008
  b = 1 / tau_c
  a = 1 - b

  xmean = randn(n,1)            # Current best mean value, should be sampled from search space...

  sigma = 1.0                   # Mutation strength, sigma (self-adapted) 

  C = covarMatrixSampler(n)

  # Calc utilities for the mu best individuals.
  utilities = utilitiesFunc(mu, lambda)

  # Now lets optimize!
  for(step in 1:num_reps)
    #################################################################
    # Decompose the covar matrix for this round
    #################################################################
    decompose!(C)

    # Generate new population
    sigmas = sigma * exp( tau * randn(1, lambda) )  # 1*lambda
    s = multivariate_normal_sample(C, n, lambda)
    z = broadcast(*, sigmas, s)                     # n*lambda
    xs = broadcast(+, xmean, z)                       # n*lambda

    # Evaluate fitness
    fitnesses = eval_fitnesses(func, xs, lambda)

    # Assign weights to the best individuals according to the utilities vector.
    weights = assign_weights(lambda, fitnesses, utilities)

    # Calculate new mean value based on weighted version of steps in population.
    xmean += (z * weights)

    # Print some info
    print("$(step): Mean point fitness = $(func(xmean)) for\n  x = $(xmean')")

    # Update the covariance matrix
    ws = broadcast(*, weights', s)
    update_covariance_matrix!(C, (ws * ws'), a)

    # Adapt sigma for next round
    sigma = sigmas * weights
  end

  return xmean, func(xmean)
end

function eval_fitnesses(func, xs, lambda = size(xs, 2))
  fitnesses = zeros(lambda)
  for(i in 1:lambda)
    fitnesses[i] = func(xs[:,i])
  end
  fitnesses
end

# We create different types of Covariance matrix samplers based on different
# decompositions.
abstract CovarianceMatrixSampler

type EigenCovarSampler <: CovarianceMatrixSampler
  C::Array{Float64,2}
  B::Array{Float64,2}
  diagD::Array{Float64,1}

  EigenCovarSampler(n) = begin
    new(eye(n,n), eye(n,n), ones(n))
  end
end

function update_covariance_matrix!(cms::CovarianceMatrixSampler, delta, a)
  C = a * cms.C + (1 - a) * delta
  cms.C = triu(C) + triu(C,1)' # Ensure C is symmetric. Should not be needed, investigate...
end

function decompose!(cms::EigenCovarSampler)
  EV, B = eig(cms.C)
  cms.B = B
  cms.diagD = sqrt(EV)
end

function multivariate_normal_sample(cms::CovarianceMatrixSampler, n, m)
  cms.B * (cms.diagD .* randn(n, m))
end

type CholeskyCovarSampler <: CovarianceMatrixSampler
  C::Array{Float64,2}
  sqrtC::Array{Float64,2}

  CholeskyCovarSampler(n) = begin
    new(eye(n,n), eye(n,n))
  end
end

function decompose!(cms::CholeskyCovarSampler)
  cms.sqrtC = chol(cms.C)'
end

function multivariate_normal_sample(cms::CholeskyCovarSampler, n, m)
  cms.sqrtC * randn(n, m)
end


function assign_weights(n, fitnesses, utilities; minimize = true)
  us_ordered = zeros(n, 1)
  perms = sortperm(fitnesses, rev = !minimize)
  for(i in 1:n)
    us_ordered[perms[i]] = utilities[i]
  end
  us_ordered
end

function rosenbrock(x)
  n = length(x)
  return( sum( 100*( x[2:n] - x[1:(n-1)].^2 ).^2 + ( x[1:(n-1)] - 1 ).^2 ) )
end

x, f = cmsa_es(3, rosenbrock; num_reps = 1000)
println("Fitness = $(f),\n  x = $(x')")