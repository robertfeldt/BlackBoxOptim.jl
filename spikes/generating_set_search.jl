# A direction generator generates the search directions to use at each step of
# a GSS search.
abstract type DirectionGenerator end

mutable struct ConstantDirectionGen <: DirectionGenerator
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
compass_search_directions(n) = ConstantDirectionGen([Matrix{Float64}(I, n,n) -Matrix{Float64}(I, n, n)])

using BlackBoxOptim

# Generating set search as described on page 21 (405) in Kolda2003 but extended
# with random ordering or frequency adapted ordering of directions.
function generating_set_search(problem;
  max_seconds = 2*numdims(problem), max_evals_per_dim = 1e5,
  ftol = 1e-7,
  delta_tol = 1e-13,
  step_size = false, x = false, random_order = false,
  freq_adapt_order = false,
  direction_gen = compass_search_directions(numdims(problem)))

  n = numdims(problem)
  ss = search_space(problem)

  step_size = 0.50 * minimum(diameters(ss))

  @assert delta_tol > 0
  @assert step_size > delta_tol

  max_evals = max_evals_per_dim * n

  directions = [Matrix{Float64}(I, n,n) -Matrix{Float64}(I, n, n)]

  x = rand_individual(ss)
  fbest = eval1(x, problem)
  num_fevals = 1

  archive = BlackBoxOptim.TopListArchive(n, 10)
  add_candidate!(archive, fbest, x[:], num_fevals)

  fc = fbest
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

    if num_fevals > max_evals
      termination_reason = "Exceeded function eval budget"
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

    for direction in order
      candidate = x + step_size .* directions[:, direction]

      fc = eval1(candidate, problem)
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
      add_candidate!(archive, fbest, x[:], num_fevals)
      println(num_fevals, ": fbest = ", fbest)

      if fitness_is_within_ftol(problem, ftol, fbest)
        termination_reason = "Within ftol"
        break
      end
    else
      step_size = step_size / 2
    end
  end

  return x, fbest, num_fevals, termination_reason, archive
end

function min_distance(point, points)
  dist = norm(point - points[:,1])
  for i in 2:size(points, 2)
    dist = minimum([dist, norm(point - points[:, i])])
  end
  dist
end

function maximize_min_distance(n, points, num_samples = 100,
  sampler = () -> randn(n,1))
  xbest = xnew = sampler()
  dist = min_distance(xnew, points)
  for i in 1:num_samples
    xnew = sampler()
    newdist = min_distance(xnew, points)
    if newdist > dist
      xbest = xnew
      dist = newdist
    end
  end
  xbest
end

function restart_gss(f, n; max_seconds = 2*n, delta_tol = 1e-7,
  step_size = 1.0,
  x = false, diverse_starting_points = false, diameter = 1.0,
  random_order = false, freq_adapt_order = false,
  known_fmin = :unknown, ftol = 1e-8,
  direction_gen = compass_search_directions(n))

  now = start_time = time()
  total_fbest = Inf
  xb = total_fevals = termination_reason = 0
  num_runs = 0

  sampler_func = () -> diameter * randn(n, 1)

  if x == false
    # Generate the first starting point.
    starting_points = sampler_func()
  else
    starting_points = x
  end

  while(now < start_time + max_seconds)
    time_left = max_seconds - (now - start_time)

    num_runs += 1

    x, fbest, num_fevals, termination_reason = generating_set_search(f, n;
      max_seconds = time_left, delta_tol = delta_tol, random_order = random_order,
      freq_adapt_order = freq_adapt_order,
      step_size = step_size, x = starting_points[:,num_runs],
      known_fmin = known_fmin, ftol = ftol,
      direction_gen = direction_gen)

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

    if diverse_starting_points
      xnew = maximize_min_distance(n, starting_points, 1000, sampler_func)
    else
      xnew = sampler_func()
    end
    starting_points = hcat(starting_points, xnew)

    now = time()
  end

  if time() > start_time + max_seconds
    termination_reason = "Exceeded time budget (elasped = $(time()-start_time))"
  end

  return xb, total_fbest, total_fevals, termination_reason
end

# Run parameter study with given parameters.
function run_with_params(of, n)
  start_time = time()
  xb, fb, fes, tr = restart_gss(of, n; known_fmin = 0.0, random_order = false)
  end_time = time()
  elapsed = end_time - start_time

  prefix_header = "StartTime,ElapsedTime,Function,N,FitnessBest,TotalFuncEvals,TerminationReason"
  println(join([prefix_header, map((i) -> "x$(i)", 1:length(xb))], ","))
  println(join([strftime("%Y%m%d-%H%M%S", start_time), elapsed, of, n, fb, fes, tr, xb], ","))
end

#of = rastrigin
#n = 5
#@time fbs, fes = repeated_runs(
#  ((n) -> restart_gss(of, n; known_fmin = 0.0, random_order = false)),
#  n);

#@time generating_set_search(sphere, 100; known_fmin = 0.0, random_order = true)
#@time generating_set_search(rosenbrock, 16; known_fmin = 0.0)
#@time restart_gss(rosenbrock, 32; known_fmin = 0.0, random_order = false)
#
#@time restart_gss(rosenbrock, 16; known_fmin = 0.0)
#@time restart_gss(xtransform(16, rosenbrock), 16; known_fmin = 0.0)
#
#@time restart_gss(cigar, 32; known_fmin = 0.0)

#of = rosenbrock
#@time ts, fbs, fes = compare_params([
#  (of, 2, false, false, false, 1.0),
#  (of, 2, false, false, false, 2.0),
#  (of, 2, false, false, true, 1.0),
#  (of, 2, false, false, true, 2.0)
#  ],
#  ((f, n, fao, ro, div, diam) ->
#    restart_gss(xtransform(n, f), n;
#      known_fmin = 0.0, freq_adapt_order = fao, random_order = ro,
#      diverse_starting_points = div, diameter = diam)),
#  5
#);
#
#
#of = griewank
#n = 2
#@time fbs, fes = repeated_runs(
#  ((n) -> restart_gss(of, n; known_fmin = 0.0, random_order = false)),
#  n);


#of = rosenbrock
#@time ts, fbs, fes = compare_params([
#  (of, 10, true, false),
#  #(of, 10, false, true),
#  #(of, 10, false, false),
#  #(of, 30, true, false),
#  #(of, 30, false, true),
#  (of, 30, false, false)
#  ],
#  ((f, n, fao, ro) ->
#    restart_gss(xtransform(n, f), n;
#      known_fmin = 0.0, freq_adapt_order = fao, random_order = ro)),
#  2
#)
#
#f = schwefel2_21
#n = 16
#
#of = deceptive_cuccu2011(30, 2)
#of = rastrigin
#n = 5
#@time fbs, fes = repeated_runs(
#  ((n) -> restart_gss(of, n; known_fmin = 0.0, random_order = false)),
#  n);
#
#
## Comparing to Cuccu2011
#of = deceptive_cuccu2011(15, 2)
#of = rastrigin
#@time ts, fbs, fes = compare_params([
#  (of, 2, true, false),
#  #(of, 2, false, true),
#  (of, 2, false, false),
#  (of, 5, true, false),
#  #(of, 5, false, true),
#  (of, 5, false, false),
#  (of, 10, true, false),
#  #(of, 10, false, true),
#  (of, 10, false, false),
#  (of, 20, true, false),
#  #(of, 20, false, true),
#  (of, 20, false, false)
#  ],
#  ((f, n, fao, ro) ->
#    restart_gss(xtransform(n, f), n;
#      known_fmin = 0.0, freq_adapt_order = fao, random_order = ro)),
#  25
#);
#
## Cuccu tells the median number of generations but not clear how many func evals
## this corresponds to. Lets assume he uses the standard 4+3*log(d) rule to set
## population size => for d=20 we should multiply by 12.
## Table IV on page 5 of Cuccu2011 should thus read:
##   d=2,  2112
##   d=5,  1872
##   d=10, 633
##   d=20, 16488
##
## These results for d=20 are impressive.
