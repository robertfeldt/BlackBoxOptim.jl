function solis_wets(problem; 
  max_evals_per_dim = 1e5, ftol = 1e-7,
  max_success_steps = 5, max_fail_steps = 3,
  rnggen = (sigma) -> Normal(0, sigma))

  start_time = time()

  n = numdims(problem)
  ss = search_space(problem)
  num_successes = num_failures = 0
  sigma = sqrt(0.30 * minimum(diameters(ss)))
  diff = bias = zeros(n)
  rng = rnggen(sigma)
  max_evals = max_evals_per_dim * n

  xbest = rand_individual(ss)
  termination_reason = "unknown"

  a = BlackBoxOptim.TopListArchive(n, 10)
  fbest = eval1(xbest, problem)
  num_fevals = 1
  add_candidate!(a, fbest, xbest[:], num_fevals)

  while(true)
    # Terminate if fitness within ftol of known fmin
    if fitness_is_within_ftol(problem, ftol, fbest)
      termination_reason = "Within ftol"
      break
    end

    if num_fevals > max_evals
      termination_reason = "Max evals budget exceeded"
      break
    end

    diff = rand(rng, n)
    xplus = xbest + bias + diff
    fxplus = eval1(xplus, problem)
    num_fevals += 1

    if fxplus < fbest

      fbest = fxplus
      add_candidate!(a, fxplus, xplus[:], num_fevals)
      num_successes += 1
      num_failures = 0
      bias = 0.2 * bias + 0.4 * (diff + bias)
      xbest = xplus
      println( "$(num_fevals): New best, fitness = $(fbest), width = $(width_of_confidence_interval(a, 0.01))")

    else

      xminus = xbest - bias - diff
      fxminus = eval1(xminus, problem)
      num_fevals += 1

      if fxminus < fbest

        fbest = fxminus
        add_candidate!(a, fxminus, xminus[:], num_fevals)
        num_successes += 1
        num_failures = 0
        bias = bias - 0.4 * (diff + bias)
        xbest = xminus
        println( "$(num_fevals): New best, fitness = $(fbest), width = $(width_of_confidence_interval(a, 0.01))")

      else
        bias = bias / 2
        num_failures += 1
        num_successes = 0
      end

    end

    if num_successes > max_success_steps
      sigma = sigma * 2
      rng = rnggen(sigma)
      num_successes = 0
    elseif num_failures > max_fail_steps
      sigma = max(1e-20, sigma / 2)
      rng = rnggen(sigma)
      num_failures = 0
    end

  end

  return xbest, best_fitness(a), num_fevals, termination_reason, a
end

#type SolisWetsOptimizer <: Optimizer
#  # Current population: 1st column is current solution, 2nd columnn is current
#  # candidate to consider. 
#  pop::Array{Float64, 2}
#  N::Int      
#  bias::Array{Float64, 1}   # Bias
#  diff::Array{Float64, 1}   # Latest Diff vector.
#  sigma::Float64            # Measure of spread in random sampling. Note! Sigma used in stead of rho=sigma^2
#  rng::ContinuousMultivariateDistribution
#  ss::SearchSpace
#  num_successes::Int
#  num_failures::Int
#  state_is_plus::Bool       # Indicates which state we are in, i.e. if next candidate to generate is the plus or the minus candidate.
#
#  SolisWetsOptimizer(ss::SearchSpace, x = false; bias = false, sigma = false, rng = false) = begin
#    if rng == false
#      sigma = sigma ? sigma : (0.5*mean(diameters(ss)))
#      rng = Normal(0, sigma)
#    end
#
#    N = numdims(ss)
#
#    pop = zeros(N, 2)
#    if x != false
#      pop[:, 1] = x
#    else
#      pop[:, 1] = rand(rng, N)
#    end
#
#    bias = zeros(N)
#    diff = zeros(N)
#
#    new(pop, N, bias, bdiff, sigma, rng, ss, 0, 0, true)
#  end
#end
#
## Get a new set of candidate solutions.
#function ask(sw::SolisWetsOptimizer)
#  if sw.state_is_plus
#    sw.diff = rand(sw.rng, sw.N)
#    xdiff = sw.bias + sw.diff
#  else
#    xdiff = -sw.bias - sw.diff
#  end
#
#  sw.pop[:, 2] = sw.pop[:, 1] + xdiff
#
#  return sw.pop, [2], [1]
#end
#
#function increase_sigma(sw::SolisWetsOptimizer)
#  sw.sigma = sw.sigma * 2
#  sw.rng = Normal(0, sw.sigma)
#end
#
#function decrease_sigma(sw::SolisWetsOptimizer)
#  sw.sigma = sw.sigma / 2
#  sw.rng = Normal(0, sw.sigma)
#end
#
#function tell(sw::SolisWetsOptimizer, pop, ranks)
#  success = ranks[1] == 2  # Success if the 2nd individual is ranked first.
#
#  if success
#    sw.pop[:, 1] = sw.pop[:, 2]
#    sw.num_successes += 1
#    sw.num_failures = 0
#    sign = sw.state_is_plus ? (+1) : (-1)
#    sw.bias = 0.2 * sw.bias + sign * 0.4 * (sw.diff + sw.bias)
#    # With a success we always go next in the plus direction.
#    sw.state_is_plus = true
#  else
#    if sw.state_is_plus == false
#      sw.num_successes = 0
#      sw.num_failures += 1
#      sw.state_is_plus = true
#    end
#  end
#
#  if sw.num_successes > 5
#    increase_sigma(sw)
#    sw.num_successes = 0
#  elseif sw.num_failures > 3
#    decrease_sigma(sw)
#    sw.num_failures = 0
#  end
#end