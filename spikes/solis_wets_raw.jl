using Distributions
using BlackBoxOptim

function solis_wets(problem; max_evals_per_dim = 1e4, ftol = 1e-8,
  max_success_steps = 5, max_fail_steps = 3,
  min_improvement_potential = 1e-8, rnggen = (sigma) -> Normal(0, sigma))
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
  fbest = eval1(xbest, p)
  add_candidate!(a, fbest, xbest[:])
  num_fevals = 1

  # We map each barrier breach log10 value to its (num_fevals, fitness, width, fip)
  barrier_breaches = Dict{Int64, (Int64, Float64, Float64, Float64)}()

  while(true)
    # Terminate if fitness within ftol of known fmin
    if fitness_is_within_ftol(problem, ftol, fbest)
      termination_reason = "Within ftol"
      break
    end

    # Save info about the fitness improvement potential each time a new log10
    # ftol barrier is breached.
    barrier_class = int(ceil(log10(fbest)))
    if !haskey(barrier_breaches, barrier_class)
      w = width_of_confidence_interval(a, 0.01)
      fip = fitness_improvement_potential(a, 0.01)
      barrier_breaches[barrier_class] = (num_fevals, fbest, w, fip)
    end

    # Terminate if improvement potential is less than 0.01% at p-value 0.01
    #if fitness_improvement_potential(a, 0.01) < min_improvement_potential
    #  termination_reason = "Not enough improvement potential"
    #  break
    #end

    if num_fevals > max_evals
      termination_reason = "Max evals budget exceeded"
      break
    end

    diff = rand(rng, n)
    xplus = xbest + bias + diff
    fxplus = eval1(xplus, p)
    num_fevals += 1

    if fxplus < fbest

      fbest = fxplus
      add_candidate!(a, fxplus, xplus[:])
      num_successes += 1
      num_failures = 0
      bias = 0.2 * bias + 0.4 * (diff + bias)
      xbest = xplus
      println( "$(num_fevals): New best, fitness = $(fbest), fip = $(fitness_improvement_potential(a, 0.01)), width = $(width_of_confidence_interval(a, 0.01))")

    else

      xminus = xbest - bias - diff
      fxminus = eval1(xminus, p)
      num_fevals += 1

      if fxminus < fbest

        fbest = fxminus
        add_candidate!(a, fxminus, xminus[:])
        num_successes += 1
        num_failures = 0
        bias = bias - 0.4 * (diff + bias)
        xbest = xminus
        println( "$(num_fevals): New best, fitness = $(fbest), fip = $(fitness_improvement_potential(a, 0.01)), width = $(width_of_confidence_interval(a, 0.01))")

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
      sigma = max(1e-7, sigma / 2)
      rng = rnggen(sigma)
      num_failures = 0
    end

  end

  return xbest, best_fitness(a), num_fevals, termination_reason, a, barrier_breaches
end

gaussianrng = (sigma) -> Normal(0, sigma)
levyrng = (sigma) -> Levy(0, sigma)

p = as_fixed_dim_problem(BlackBoxOptim.example_problems["Rosenbrock"], 2^6)
problem = BlackBoxOptim.shifted(p)
xb, fb, nf, r, a, bbs = solis_wets(problem; min_improvement_potential = 1e-8, 
  max_fail_steps = 3, max_evals_per_dim = 1e6, rnggen = gaussianrng)

println("reason = $(r), fitness = $(fb)")
bkeys = sort(collect(keys(bbs)), rev = true)
for(k in bkeys)
  println("$(k): $(bbs[k])")
end