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
function cmsa_es(problem; 
  max_seconds = 2*numdims(problem), max_evals_per_dim = 1e5,
  ftol = 1e-7,
  max_rounds_without_improvement = 1000,
  mu = false, lambda = false, 
  covarMatrixSampler = EigenCovarSampler,
  utilitiesFunc = linear_utilities,

  trace = true)

  n = numdims(problem)
  ss = search_space(problem)

  max_evals = max_evals_per_dim * n

  if lambda == false
    # In the paper they try population sizes 8, 4n and 4*n*n, lets try
    # in between:
    #lambda = int(8 + (4*n*n - 8) * rand())
    lambda = 4*n
  end

  if mu == false
    mu = int(maximum([ceil(lambda/2), 2.0]))
  end

  tau = 1 / sqrt(2*n)           # Equation (1) on page 5 in Beyer2008 
  tau_c = 1 + n * (n + 1) / 2   # Equation (2) on page 5 in Beyer2008
  a = 1 - 1 / tau_c

  sigma = 1.0                   # Mutation strength, sigma (self-adapted) 
  C = covarMatrixSampler(n)
  utilities = utilitiesFunc(mu, lambda)

  xbest = xmean = rand_individual(ss)   # Current best mean value.
  fbest = eval1(xbest, problem)
  num_fevals = 1

  archive = BlackBoxOptim.TopListArchive(n, 10)
  add_candidate!(archive, fbest, xbest[:], num_fevals)

  fevals_last_best = 0
  next_print_covar = 100
  termination_reason = "?"

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
    fitnesses = eval_fitnesses(problem, xs, lambda)
    num_fevals += lambda

    # Check if best new fitness is best ever and print some info if tracing.
    indbest = indmin(fitnesses)
    fbest_new = fitnesses[indbest]
    add_candidate!(archive, fbest_new, xs[:, indbest], num_fevals)

    if fbest_new < fbest
      xbest = xs[:, indbest]
      fbest = fbest_new
      fevals_last_best = num_fevals

      if fitness_is_within_ftol(problem, ftol, fbest)
        termination_reason = "Within ftol"
        break
      end

      if trace
        println("$(num_fevals): Best fitness = $(fbest)")
        if num_fevals > next_print_covar
          next_print_covar = num_fevals + 100
          #println("covar matrix:")
          #show(C.C)
          println("\ncovar summary: ", sumstats(C.C, (x) -> @sprintf("%.2e", x)))
        end
      end
    end

    # Assign weights to the best individuals according to the utilities vector.
    weights = assign_weights(lambda, fitnesses, utilities)

    # Calculate new mean value based on weighted version of steps in population.
    xmean += (z * weights)

    # Update the covariance matrix
    ws = broadcast(*, weights', s)
    update_covariance_matrix!(C, (ws * ws'), a)

    # Terminate if max covar value is very small
    if maximum(C.C) < 1e-300
      termination_reason = "Max covar value very small"
      break
    end

    # Terminate if too long since no improvement
    if (num_fevals - fevals_last_best) > max_rounds_without_improvement*lambda
      termination_reason = "Too long since improvement"
      break
    end

    # Adapt sigma for next round
    sigma = sigmas * weights
  end

  # Save fitness history to a csv file
  csvfile = strftime("cmsa_es_%Y%m%d_%H%M%S.csv", start_time)
  save_fitness_history_to_csv_file(archive, csvfile)
  println("Saved fitness history to file: $(csvfile)")

  return xbest, fbest, num_fevals, termination_reason, archive
end

function eval_fitnesses(problem, xs, lambda = size(xs, 2))
  fitnesses = zeros(lambda)
  for(i in 1:lambda)
    fitnesses[i] = eval1(xs[:,i], problem)
  end
  fitnesses
end

function assign_weights(n, fitnesses, utilities; minimize = true)
  us_ordered = zeros(n, 1)
  perms = sortperm(fitnesses, rev = !minimize)
  for(i in 1:n)
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
