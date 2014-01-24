# This is Subset Tournament CMSA-ES as proposed by Robert Feldt in the paper:
#  R. Feldt, "Covariate Subset Tournaments for High-Dimensional Blackbox Optimization with Covariance Matrix Adapting Evolutionary Strategies", 2014

function normalize_utilities(utilities)
  utilities / sum(utilities)
end

function linear_utilities(mu, lambda)
  normalize_utilities(hcat(ones(1, mu), zeros(1, lambda-mu)))
end

function log_utilities(mu, lambda)
  normalize_utilities(hcat((log(mu+1) - log(1:mu))', zeros(1, lambda-mu)))
end

# Optimize the n-dimensional objective func with a (mu,lambda) CMSA-ES using
# subset tournaments to optimize subsets of variables in rounds.
function st_cmsa_es(p; 
  trace = true,

  # Stopping criteria related
  max_seconds = 4*numdims(p), max_evals_per_dim = 1e7,
  ftol = 1e-7, xtol = 1e-10, stol = 1e-10,
  max_rounds_without_improvement = 200,

  # Starting points, will be random unless specified
  xmean = false,

  # Algorithm specific params:
  covarMatrixSampler = SubsetCholeskyCovarSampler,
  utilitiesFunc = log_utilities,
  lambda = 4*numdims(p),
  mu = int(max(ceil(lambda/rand(4:20)), 1.0)),
  tau = 1 / sqrt(2*numdims(p)), # Equation (1) on page 5 in Beyer2008
  tau_c = 1 + numdims(p) * (numdims(p) + 1) / (2 * mu),  # Equation (2) on page 5 in Beyer2008
  sigma = 0.05*rand(1:8)*minimum(diameters(search_space(p))),
  decompose_covar_prob = 0.4,

  # Subsetting specific parameters:
  subset_size = int(floor(0.20*numdims(p))),
  subsets_per_tournament = 2,
  num_rounds_per_tournament = 1,
  num_rounds_of_optimization_between_tournaments = 30,
  subset_selection_mechanism = :random

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

  # Init for subset tournaments
  st_state = :optimization   # Can be either :tournament or :optimization and starts with :optimization to ensure we set up a new tournament
  num_optimization_rounds_for_this_subset = num_rounds_of_optimization_between_tournaments # Also ensures we are end of opt round => new tournament will be set up
  st_current_subset = 1  # When in :tournament mode this is the index to the subset currently being evaluated, can be in range 1-subsets_per_tournament
  st_subsets = Array{Int64, 1}[] # Subsets that are currently in a tournament
  num_tournament_rounds_for_this_subset = 0   # Number of eval rounds we have ran with current subset in tournament mode
  fitness_per_tournament_round = zeros(num_rounds_per_tournament, subsets_per_tournament)

  # Keep stats per covariate being optimized for how effective the optimization is when they are included.
  # We save the expected relative change per round of optimization and use a history parameter in [0.0, 1.0]
  # to decide how much history should be saved when updating it, i.e. newval = oldvalue * h + delta * (1-h).
  # These stats values are used in multiple tournaments like so:
  #   1. Repeat until the right number of vars has been added to selected set
  #   2.   Randomly sample two unselected vars
  #   3.   Select the one with the best stats value (lowest value if minimizing)
  improvement_per_var = 0.01 * ones(N) # Start with some value since it will be adapted as we go...

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

    if st_state == :tournament
      if num_tournament_rounds_for_this_subset < num_rounds_per_tournament
        num_tournament_rounds_for_this_subset += 1
      else
        st_current_subset += 1
        if st_current_subset <= subsets_per_tournament
          set_subset!(C, st_subsets[st_current_subset])
          num_tournament_rounds_for_this_subset = 1
        else # No more tournaments needed, select best subset and start optimizing with it
          winning_subset = select_winning_subset(st_subsets, fitness_per_tournament_round)
          set_subset!(C, st_subsets[winning_subset])
          st_current_subset = winning_subset
          st_state = :optimization
          num_optimization_rounds_for_this_subset = 1
        end
      end
    else # In optimization mode
      if num_optimization_rounds_for_this_subset < num_rounds_of_optimization_between_tournaments
        num_optimization_rounds_for_this_subset += 1
      else
        # We have now used the current subset for a full set of optimization rounds so we set up for a new tournament
        st_subsets = select_new_subsets(st_subsets, N, subsets_per_tournament, st_current_subset, subset_selection_mechanism, subset_size)
        st_state = :tournament
        st_current_subset = 1
        set_subset!(C, st_subsets[st_current_subset])
        num_tournament_rounds_for_this_subset = 1
      end
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

    # Save info about the fitnesses if we are in tournament mode
    if st_state == :tournament
      fitness_per_tournament_round[num_tournament_rounds_for_this_subset, st_current_subset] = fbest_new
    end

    if fbest_new < fbest
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
          #println("covar summary: ", sumstats(C.C, (x) -> @sprintf("%.2e", x)))
          println("sigma: ", sigma)
        end
      end
    else
      if (num_fevals - fevals_last_best) > max_rounds_without_improvement * lambda
        termination_reason = "Max rounds without improvement reached"
        break
      end
    end

    # Assign weights to the best individuals according to the utilities vector.
    weights = assign_weights(lambda, fitnesses, utilities)

    # Calculate new mean value based on weighted version of steps in population.
    xmean += (z * weights)

    # Update the covariance matrix
    uc = zeros(Float64, N, N)
    for(i in 1:lambda)
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

    # Terminate if sigma very small
    #println("sigma = $(sigma[1,1])")
    if sigma[1,1] < stol
      termination_reason = "Sigma too small"
      break
    end
  end

  return xbest, fbest, num_fevals, termination_reason, archive
end

function select_winning_subset(subsets, fitness_per_tournament_round)
  if length(size(fitness_per_tournament_round)) > 1
    fitness_summary = mean(fitness_per_tournament_round, 1)
    indmin(fitness_summary)
  else
    indmin(fitness_per_tournament_round)
  end
end

function select_new_subsets(subsets, n, num_subsets, winning_subset, subset_selection_mechanism, subset_size)
  if length(subsets) == 0 || subset_selection_mechanism == :random

    generate_random_subsets(num_subsets, subset_size, n)

  else # We should keep the winner

    # Generate new subsets randomly
    new_subsets = generate_random_subsets(num_subsets-1, subset_size, n)

    # Add the previous winning subset
    push!(new_subsets, subsets[winning_subset])
    new_subsets

  end
end

# Select subsets based on running binary tournaments between stats values
function select_new_subsets_based_on_stats(n, num_subsets, subset_size, stats, smaller_stats_is_better)
  [generate_subset_based_on_stats(subset_size, n, stats, smaller_stats_is_better) for i in 1:num_subsets]
end

function generate_subset_based_on_stats(subset_size, n, stats, smaller_stats_is_better = false)
  subset = Int64[]
  candidates = shuffle(collect(1:n))
  len = n-1
  op = smaller_stats_is_better ? < : >

  for i in 1:subset_size
    index = rand(1:len)
    candidate1 = candidates[index]
    candidate2 = candidates[index+1]
    println("Comparing $(candidate1) and $(candidate2)")
    if op( stats[candidate1], stats[candidate2] )
      splice!(candidates, index)
      push!(subset, candidate1)
    else
      splice!(candidates, index+1)
      push!(subset, candidate2)
    end
    len -= 1
  end

  subset
end

function generate_random_subsets(num_subsets, subset_size, n)
  new_subsets = Array{Int64, 1}[]
  for i in 1:num_subsets
    push!(new_subsets, sort(shuffle(collect(1:n))[1:subset_size]))
  end
  new_subsets
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

# The subset CholeskyCovarSampler only uses a subset of the variables in the 
# (expensive) cholesky decomposition, the other variables are kept constant.
# However, the covariance matrix itself is always updated and saved in full
# so that the overall learning of the shape of the fitness landscape is not lost.
type SubsetCholeskyCovarSampler <: CovarianceMatrixSampler
  C::Array{Float64,2}
  sqrtC::Array{Float64,2}
  subset::Array{Int64, 1} # Indices of the currently active subset of variables

  SubsetCholeskyCovarSampler(n) = begin
    new(eye(n,n), eye(n,n), collect(1:n))
  end
end

# Set new subset.
function set_subset!(cms::SubsetCholeskyCovarSampler, subset)
  cms.subset = subset
  decompose!(cms)
end

using JSON

function decompose!(cms::SubsetCholeskyCovarSampler)
  try 
    subset_C = cms.C[cms.subset, cms.subset]
    cms.sqrtC = chol(subset_C)'
  catch error
    # We don't update if there is some problem
    println("ERROR: Could not do cholesky decomposition!")
    show(error)
  end
end

function multivariate_normal_sample(cms::SubsetCholeskyCovarSampler, n, m)
  subset_size = length(cms.subset)

  subset_inv_c = 1

  # Calc the inverted covar matrix only for the subset
  try 
    subset_inv_c = cms.sqrtC * randn(subset_size, m)
  catch error
    println("Error in subset_inv_c multiplication, size(cms.sqrtC) = $(size(cms.sqrtC)), subset_size = $(subset_size), m = $(m), size(cms.C) = $(size(cms.C))")
    show(error)
  end

  # The rest should be zero. Maybe faster if we do this with a sparse matrix
  # since most of them will be zero?
  inv_c = zeros(n, m)
  for i in 1:subset_size
    si = cms.subset[i]
    for j in 1:m
      inv_c[si, j] = subset_inv_c[i, j]
    end
  end

  inv_c
end
