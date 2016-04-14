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
function cmsa_es(p;
  trace = true,

  # Stopping criteria related
  max_seconds = 4*numdims(p), max_evals_per_dim = 1e7,
  ftol = 1e-7, xtol = 1e-10, stol = 1e-10,
  max_rounds_without_improvement = 500,

  # Starting points, will be random unless specified
  xmean = false,

  # Algorithm specific params:
  covarMatrixSampler = CholeskyCovarSampler,
  utilitiesFunc = log_utilities,
  lambda = 4*numdims(p),
  mu = int(max(ceil(lambda/rand(4:20)), 1.0)),
  tau = 1 / sqrt(2*numdims(p)), # Equation (1) on page 5 in Beyer2008
  tau_c = 1 + numdims(p) * (numdims(p) + 1) / (2 * mu),  # Equation (2) on page 5 in Beyer2008
  sigma = 0.05*rand(1:8)*minimum(diameters(search_space(p))),
  decompose_covar_prob = 0.4,
  max_successes_before_increasing = 3,
  max_failures_before_decreasing = 3

  )

  N = numdims(p)
  ss = search_space(p)
  max_evals = max_evals_per_dim * N
  a = 1 - 1 / tau_c
  C = covarMatrixSampler(N)
  utilities = utilitiesFunc(mu, lambda)
  xbest = xmean = rand_individual(ss)   # Current best mean value.
  fbest = eval1(xbest, p)
  num_fevals = 1
  archive = BlackBoxOptim.TopListArchive(N, 10)
  add_candidate!(archive, fbest, xbest[:], num_fevals)
  fevals_last_best = num_fevals
  next_print_covar = 100
  termination_reason = "?"

  num_successes = 0
  num_failures = 0

  start_time = time()

  # Now lets optimize! Ensure we run at least one iteration.
  while(true)

    if (time() - start_time) > max_seconds
      termination_reason = "Exceeded time budget"
      break
    end

    if num_fevals > max_evals
      termination_reason = "Exceeded function eval budget"
      break
    end

    # Decompose only with a given probability => saves time
    if rand() <= decompose_covar_prob
      decompose!(C)
    end

    # Generate new population
    sigmas = sigma * exp( tau * randn(1, lambda) )  # 1*lambda
    s = multivariate_normal_sample(C, N, lambda)
    z = broadcast(*, sigmas, s)                     # n*lambda
    xs = repmat(xmean, 1, lambda) + z               # n*lambda

    # Evaluate fitness
    fitnesses = eval_fitnesses(p, xs, lambda)
    num_fevals += lambda

    # Check if best new fitness is best ever and print some info if tracing.
    indbest = indmin(fitnesses)
    fbest_new = fitnesses[indbest]

    if fbest_new < fbest
      num_successes += 1 # From Solis-Wets
      num_failures = 0   # From Solis-Wets
      xbest = xs[:, indbest]
      fbest = fbest_new
      add_candidate!(archive, fbest, xbest, num_fevals)
      fevals_last_best = num_fevals

      if fitness_is_within_ftol(p, ftol, fbest)
        termination_reason = "Within ftol"
        break
      end

      if trace
        println("$(num_fevals): Best fitness = $(fbest)")
        if num_fevals > next_print_covar
          next_print_covar = num_fevals + 100
          println("covar summary: ", sumstats(C.C, (x) -> @sprintf("%.2e", x)))
          println("sigma: ", sigma)
        end
      end
    else
      if (num_fevals - fevals_last_best) > max_rounds_without_improvement * lambda
        termination_reason = "Max rounds without improvement reached"
        break
      end
      num_failures += 1 # From Solis-Wets
      num_successes = 0   # From Solis-Wets
    end

    # Assign weights to the best individuals according to the utilities vector.
    weights = assign_weights(lambda, fitnesses, utilities)

    # Calculate new mean value based on weighted version of steps in population.
    xmean += (z * weights)

    # Update the covariance matrix
    uc = zeros(Float64, N, N)
    for i in 1:lambda
      if weights[i] > 0.0
        se = s[:,i]
        uc += weights[i] * (se * se')
      end
    end
    update_covariance_matrix!(C, uc, a)
    #ws = broadcast(*, weights', s)
    #update_covariance_matrix!(C, (ws * ws'), covar_learning_rate)

    # Adapt sigma for next round
    sigma = sigmas * weights

    # From solis-wets
    if (max_successes_before_increasing >= 1.0) && (num_successes >= max_successes_before_increasing)
      sigma = sigma * 2
      num_successes = 0
    elseif (max_failures_before_decreasing >= 1.0) && (num_failures >= max_failures_before_decreasing)
      sigma = sigma / 2
      num_failures = 0
    end

    # Terminate if sigma very small
    #println("sigma = $(sigma[1,1])")
    if sigma[1,1] < stol
      termination_reason = "Sigma too small"
      break
    end
  end

  return xbest, fbest, num_fevals, termination_reason, archive
end

function restart_cmsa_es(p;
  trace = true,

  # Stopping criteria related
  max_seconds = 2*numdims(p), max_evals_per_dim = 1e7,
  ftol = 1e-7, xtol = 1e-10, stol = 1e-10,

  # Starting points, will be random unless specified
  xmean = false,

  # Start values for algorithm specific params:
  lambda = 4*numdims(p),
  lambda_divisor = 7,
  covarMatrixSampler = CholeskyCovarSampler,
  utilitiesFunc = log_utilities,
  tau = 1 / sqrt(2*numdims(p)), # Equation (1) on page 5 in Beyer2008
  covar_learning_rate = 0.88,

  )

  now = start_time = time()

  total_fbest = Inf
  xb = total_fevals = termination_reason = 0
  num_runs = 0

  while(now < start_time + max_seconds)
    time_left = max_seconds - (now - start_time)

    num_runs += 1

    x, fbest, num_fevals, termination_reason, archive = cmsa_es(p;
      trace = trace, max_seconds = time_left,
      ftol = ftol, xtol = xtol, stol = stol,
      max_evals_per_dim = max_evals_per_dim,
      max_rounds_without_improvement = 1000,
      lambda = lambda,
      mu = int(maximum([ceil(lambda/lambda_divisor), 2.0])),
      covarMatrixSampler = CholeskyCovarSampler, utilitiesFunc = log_utilities,
      tau = tau, covar_learning_rate = covar_learning_rate
    );

    lambda += (4*numdims(p))
    lambda_divisor = 5 + rand(1:15)
    sigma = 0.10*rand(1:4)*minimum(diameters(search_space(p)))
    covar_learning_rate = 0.85+0.01*rand(1:4)

    total_fevals += num_fevals

    if fbest < total_fbest
      total_fbest = fbest
      xb = x
      println(num_runs, ": New global best = $(total_fbest)")

      if fitness_is_within_ftol(p, ftol, fbest)
        termination_reason = "Within ftol"
        break
      end
    end

    now = time()
  end

  if time() > start_time + max_seconds
    termination_reason = "Exceeded time budget (elasped = $(time()-start_time))"
  end

  return xb, total_fbest, total_fevals, termination_reason, archive
end

function eval_fitnesses(problem, xs, lambda = size(xs, 2))
  fitnesses = zeros(lambda)
  for i in 1:lambda
    fitnesses[i] = eval1(xs[:,i], problem)
  end
  fitnesses
end

function assign_weights(n, fitnesses, utilities; minimize = true)
  us_ordered = zeros(n, 1)
  perms = sortperm(fitnesses, rev = !minimize)
  for i in 1:n
    us_ordered[perms[i]] = utilities[i]
  end
  us_ordered
end

# We create different types of Covariance matrix samplers based on different
# decompositions.
abstract CovarianceMatrixSampler

function update_covariance_matrix!(cms::CovarianceMatrixSampler, delta, a)
  C = a * cms.C + (1 - a) * delta
  cms.C = triu(C) + triu(C,1)' # Ensure C is symmetric. Should not be needed, investigate...
end

type EigenCovarSampler <: CovarianceMatrixSampler
  C::Array{Float64,2}
  B::Array{Float64,2}
  diagD::Array{Float64,1}

  EigenCovarSampler(n) = begin
    new(eye(n,n), eye(n,n), ones(n))
  end
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

# Block-wise Covar sampling selects a block of coordinates for which the
# matrix is decomposed but then discards the covariance for the coordinates
# not in the selected block and only uses their (individual) variances.
# Since this is simpler to relate to the Eigen decomposition we start with
# that version.
