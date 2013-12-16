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

# This weights both positive and negative ones
function active_linear_utilities(mu, lambda)
  pos = ones(1, mu) / mu
  rest = lambda-mu
  neg = (-1) * ones(1, rest) / (2*rest)
  2 * hcat(pos, neg)'
end

# Optimize the n-dimensional objective func with a (mu,lambda) CMSA-ES.
function cmsa_es(n, func; 
  mu = false, lambda = false, max_fevals = 2000,
  covarMatrixSampler = EigenCovarSampler,
  utilitiesFunc = linear_utilities,
  trace = true)

  if lambda == false
    # In the paper they try population sizes 8, 4n and 4*n*n, lets try
    # in between:
    #lambda = int(8 + (4*n*n - 8) * rand())
    lambda = 4*n
    # We should probably use an IPOP scheme for restarts...
  end

  if mu == false
    mu = max(int(log(lambda)), 2)
  end

  tau = 1 / sqrt(2*n)           # Equation (1) on page 5 in Beyer2008 
  tau_c = 1 + n * (n + 1) / 2   # Equation (2) on page 5 in Beyer2008
  a = 1 - 1 / tau_c
  xmean = randn(n,1)            # Current best mean value, should be sampled from search space...
  sigma = 1.0                   # Mutation strength, sigma (self-adapted) 
  C = covarMatrixSampler(n)
  utilities = utilitiesFunc(mu, lambda)

  num_fevals = 0

  # Now lets optimize! Ensure we run at least one iteration.
  while(num_fevals < max_fevals)
    #################################################################
    # Decompose the covar matrix for this round
    #################################################################
    decompose!(C)

    # Generate new population
    sigmas = sigma * exp( tau * randn(1, lambda) )  # 1*lambda
    s = multivariate_normal_sample(C, n, lambda)
    z = broadcast(*, sigmas, s)                     # n*lambda
    xs = repmat(xmean, 1, lambda) + z               # n*lambda

    # Evaluate fitness
    fitnesses = eval_fitnesses(func, xs, lambda)
    num_fevals += lambda

    # Assign weights to the best individuals according to the utilities vector.
    weights = assign_weights(lambda, fitnesses, utilities)

    # Calculate new mean value based on weighted version of steps in population.
    xmean += (z * weights)

    if trace
      print("$(num_fevals): Best fitness = $(func(xmean))\n")
    end

    # Update the covariance matrix
    ws = broadcast(*, weights', s)
    update_covariance_matrix!(C, (ws * ws'), a)

    # Adapt sigma for next round
    sigma = sigmas * weights
  end

  return xmean, func(xmean), num_fevals
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
  try
    EV, B = eig(cms.C)
    cms.B = B
    cms.diagD = sqrt(EV)
  catch
    # We don't update if there is some problem
  end
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
  try 
    cms.sqrtC = chol(cms.C)'
  catch error
    # We don't update if there is some problem
  end
end

function multivariate_normal_sample(cms::CholeskyCovarSampler, n, m)
  cms.sqrtC * randn(n, m)
end

# This is not a general solution but it works for specifically my use case and
# when the operator is *
function Base.broadcast{Tv,Ti}(op, v::Array{Tv,2}, A::SparseMatrixCSC{Tv,Ti})
  I, J = findn(A)
  V = zeros(nnz(A))
  vn, vm = size(v)
  if vn >= 1 && vm == 1
    for(l in 1:nnz(A))
      row = I[l]
      V[l] = op(v[row], A.nzval[l])
    end
  elseif vn == 1 && vm >= 1
    for(l in 1:nnz(A))
      col = J[l]
      V[l] = op(v[col], A.nzval[l])
    end
  else
    throw(ArgumentError("invalid dimensions"))
  end
  sparse(I, J, V)
end

type SparseCholeskyCovarSampler <: CovarianceMatrixSampler
  C::SparseMatrixCSC{Float64,Int64}
  sqrtC::SparseMatrixCSC{Float64,Int64}

  SparseCholeskyCovarSampler(n) = begin
    new(speye(n,n), speye(n,n))
  end
end

function update_covariance_matrix!(cms::SparseCholeskyCovarSampler, delta, a)
  C = a * cms.C + (1 - a) * delta
  cms.C = C # triu(C) + triu(C,1)' # Ensure C is symmetric. Should not be needed, investigate...
end

function decompose!(cms::SparseCholeskyCovarSampler)
  try
    #println("droptol!(C)")
    #cms.C = Base.droptol!(cms.C, 1e-4)
    #println("nnz(C) = $(nnz(cms.C)), min = $(minimum(cms.C.nzval)), max = $(maximum(cms.C.nzval))")
    #println("cholfact(C)")
    t = cholfact(cms.C)
    #println("sparse(t)")
    s = sparse(t)
    #println("s'")
    sqrtC = s'
    #println("droptol!(sqrtC)")
    #cms.sqrtC = Base.droptol!(sqrtC, 1e-4)
    #show(cms.sqrtC)
  catch error
    # We don't update if there is some problem
    println("ERROR: ", error)
  end
end

function multivariate_normal_sample(cms::SparseCholeskyCovarSampler, n, m)
  #if n*m > 100
    # We select a random density in [0.10, 0.60]
  #  density = 0.10 + 0.60 * rand()
  #else
  #  density = 0.90
  #end
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

#function rosenbrock(x)
#  n = length(x)
#  sum( 100*( x[2:n] - x[1:(n-1)].^2 ).^2 + ( x[1:(n-1)] - 1 ).^2 )
#end
#
#function sphere(x)
#  sum(x .^ 2)
#end
#
#dim = 64
#tic()
#x, f, fevals = cmsa_es(dim, rosenbrock; max_fevals = max_fevals, covarMatrixSampler = SparseCholeskyCovarSampler, trace = true)
#x, f, fevals = cmsa_es(dim, rosenbrock; max_fevals = int(1e7), covarMatrixSampler = CholeskyCovarSampler, lambda = 4*dim*dim, trace = true)
#t = toq()
#println("fitness = $(f)")
#println("time = $(t)")