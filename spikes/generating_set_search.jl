# A direction generator generates the search directions to use at each step of 
# a GSS search.
abstract DirectionGenerator

type ConstantDirectionGen <: DirectionGenerator
  directions::Array{Float64, 2}

  ConstantDirectionGen(directions) = begin
    new(directions)
  end
end

function directions_for_k(cg::ConstantDirectionGen, k)
  cg.directions # Always the same regardless of k...
end

# We can easily do a compass search with GSS by generating directions 
# individually (+ and -) for each coordinate.
compass_search_directions(n) = ConstantDirectionGen([eye(n,n) -eye(n, n)])

using BlackBoxOptim

# Generating set search as described on page 21 (405) in Kolda2003 but extended
# with random ordering or frequency adapted ordering of directions.
function generating_set_search(f, n; max_seconds = 2*n, delta_tol = 1e-13, 
  step_size = 1.0, x = false, random_order = false,
  freq_adapt_order = false,
  known_fmin = :unknown, ftol = 1e-8,
  direction_gen = compass_search_directions(n))

  @assert delta_tol > 0
  @assert step_size > delta_tol

  if x == false
    x = step_size * randn(n, 1)
  end

  directions = [eye(n,n) -eye(n, n)]

  num_fevals = 0
  fc = fbest = Inf
  candidate = zeros(n, 1)
  k = 0
  termination_reason = :unknown
  order = collect(1:size(directions, 2))

  if freq_adapt_order
    fa_directions = BlackBoxOptim.FrequencyAdapter(size(directions, 2))
  end

  start_time = time()

  while(true)

    if step_size < delta_tol
      termination_reason = "Step size smaller than delta_tol"
      break
    end

    if (time() - start_time) > max_seconds
      termination_reason = "Exceeded time budget"
      break
    end

    # Get the directions for this iteration
    k += 1
    directions = directions_for_k(direction_gen, k)

    # Check all directions to find a better point
    found_better = false

    # Freq adaptation takes precedence over random ordering so check first.
    if freq_adapt_order
      order = BlackBoxOptim.create_new_block!(fa_directions)
    elseif random_order == true
      shuffle!(order)
    end

    for(direction in order)
      candidate = x + step_size .* directions[:, direction]

      fc = f(candidate)
      num_fevals += 1

      # Update freq based on progress
      if freq_adapt_order
        improvement = (fbest - fc) / fbest
        #BlackBoxOptim.update!(fa_directions, direction, max(0, improvement))
        BlackBoxOptim.update!(fa_directions, direction, improvement)
      end

      if fc < fbest
        found_better = true
        break
      end

    end

    if found_better
      fbest = fc
      x = candidate
      println(num_fevals, ": fbest = ", fbest)

      if known_fmin != :unknown && abs(fbest - known_fmin) < ftol
        termination_reason = "Within ftol"
        break
      end
    else
      step_size = step_size / 2
    end
  end

  return x, fbest, num_fevals, termination_reason
end

function restart_gss(f, n; max_seconds = 2*n, delta_tol = 1e-7, 
  step_size = 1.0, x = false, random_order = false,
  freq_adapt_order = false,
  known_fmin = :unknown, ftol = 1e-8,
  direction_gen = compass_search_directions(n))

  now = start_time = time()
  total_fbest = Inf
  xb = total_fevals = termination_reason = 0
  num_runs = 0

  while(now < start_time + max_seconds)
    time_left = max_seconds - (now - start_time)

    x, fbest, num_fevals, termination_reason = generating_set_search(f, n; 
      max_seconds = time_left, delta_tol = delta_tol, random_order = random_order,
      freq_adapt_order = freq_adapt_order,
      step_size = step_size, x = false, known_fmin = known_fmin, ftol = ftol,
      direction_gen = direction_gen)

    num_runs += 1
    step_size = step_size * 2
    #delta_tol /= 10

    total_fevals += num_fevals

    if fbest < total_fbest
      total_fbest = fbest
      xb = x
      println(num_runs, ": New global best = $(total_fbest)")
      if termination_reason == "Within ftol"
        break
      end
    end

    now = time()
  end

  if time() > start_time + max_seconds
    termination_reason = "Exceeded time budget (elasped = $(time()-start_time))"
  end

  return xb, total_fbest, total_fevals, termination_reason
end

# A GSS variant with an initial Levy flight search step to investigate if this
# can speed up convergence on multi-modal fitness functions.
function gss_levy()
end

@time generating_set_search(sphere, 100; known_fmin = 0.0, random_order = true)
@time generating_set_search(rosenbrock, 16; known_fmin = 0.0)
@time restart_gss(rosenbrock, 32; known_fmin = 0.0, random_order = false)

@time restart_gss(rosenbrock, 16; known_fmin = 0.0)
@time restart_gss(xtransform(16, rosenbrock), 16; known_fmin = 0.0)

@time restart_gss(cigar, 32; known_fmin = 0.0)

function sumstats(v, f = (x) -> @sprintf("%.2f", x))
  "median = $(f(median(v))), mean = $(f(mean(v))) +/- $(f(std(v)))"
end

function format_time(t)
  if t < 5e-1
    @sprintf("%.2f ms", t*1e3)
  elseif t < 30.0
    @sprintf("%.2f s", t)
  elseif t < 30*60
    @sprintf("%.2f min", t/60.0)
  else
    @sprintf("%.2f hours", t/3600.0)
  end
end

function repeated_runs(searchf, n = 16, num_runs = 10)
  fevals = zeros(num_runs)
  fbests = zeros(num_runs)
  times = zeros(num_runs)
  xs = Any[]
  for(i in 1:num_runs)
    tic()
    x, fbests[i], fevals[i], termination_reason = searchf(n)
    times[i] = toq()
    push!(xs, x)
  end

  println("\nFitness: ", sumstats(fbests))
  println("Time: ", sumstats(times, format_time))
  println("Num. evals: ", sumstats(fevals, int))
  println("")

  return fbests, fevals
end

function compare_params(params, searchf, num_runs = 10)
  num_configs = length(params)
  times = zeros(num_runs, num_configs)
  fbests = zeros(num_runs, num_configs)
  fevals = zeros(num_runs, num_configs)

  for(r in 1:num_runs)
    for(i in 1:num_configs)
      tic()
      xb, fb, fe, reason = searchf(params[i]...)
      times[r,i] = toq()
      fbests[r, i] = fb
      fevals[r, i] = float(fe)
    end
  end

  println("\n\nResults per parameter config")
  println("----------------------------")
  for(i in 1:num_configs)
    print(i, ". "); show(params[i]); println("")
    println("Fitness: ", sumstats(fbests[:,i]))
    println("Time: ", sumstats(times[:,i]))
    println("Num. evals: ", sumstats(fevals[:,i]))
    println("")
  end

  return times, fbests, fevals
end

of = rosenbrock
@time ts, fbs, fes = compare_params([
  (of, 10, true, false), 
  (of, 10, false, true), 
  (of, 10, false, false),
  (of, 30, true, false), 
  (of, 30, false, true), 
  (of, 30, false, false)
  ],
  ((f, n, fao, ro) -> 
    restart_gss(xtransform(n, f), n; 
      known_fmin = 0.0, freq_adapt_order = fao, random_order = ro)),
  3
)

f = schwefel2_21
n = 16

of = deceptive_cuccu2011(30, 2)
of = rastrigin
n = 5
@time fbs, fes = repeated_runs(
  ((n) -> restart_gss(of, n; known_fmin = 0.0, random_order = false)),
  n);


# Comparing to Cuccu2011
of = deceptive_cuccu2011(15, 2)
of = rastrigin
@time ts, fbs, fes = compare_params([
  (of, 2, true, false), 
  (of, 2, false, true), 
  (of, 2, false, false),
  (of, 5, true, false), 
  (of, 5, false, true), 
  (of, 5, false, false),
  (of, 10, true, false), 
  (of, 10, false, true), 
  (of, 10, false, false),
  (of, 20, true, false), 
  (of, 20, false, true), 
  (of, 20, false, false)
  ],
  ((f, n, fao, ro) -> 
    restart_gss(xtransform(n, f), n; 
      known_fmin = 0.0, freq_adapt_order = fao, random_order = ro)),
  25
);

# Cuccu tells the median number of generations but not clear how many func evals
# this corresponds to. Lets assume he uses the standard 4+3*log(d) rule to set
# population size => for d=20 we should multiply by 12.
# Table IV on page 5 of Cuccu2011 should thus read:
#   d=2,  2112
#   d=5,  1872
#   d=10, 633
#   d=20, 16488
#
# These results for 