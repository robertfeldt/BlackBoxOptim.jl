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

# Generating set search as described on page 21 (405) in Kolda2003:
function generating_set_search(f, n; max_seconds = 2*n, delta_tol = 1e-10, 
  step_size = 1.0, x = false, 
  known_fmin = :unknown, ftol = 1e-8,
  direction_gen = compass_search_directions(n))

  @assert delta_tol > 0
  @assert step_size > delta_tol

  if x == false
    x = randn(n, 1)
  end

  directions = [eye(n,n) -eye(n, n)]

  num_fevals = 0
  fc = fbest = Inf
  candidate = zeros(n, 1)
  k = 0
  termination_reason = :unknown

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
    for(direction in 1:size(directions, 2))
      candidate = x + step_size .* directions[:, direction]

      fc = f(candidate)
      num_fevals += 1

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
  step_size = 1.0, x = false, 
  known_fmin = :unknown, ftol = 1e-8,
  direction_gen = compass_search_directions(n))

  now = start_time = time()
  total_fbest = Inf
  xb = total_fevals = termination_reason = 0
  num_runs = 0

  while(now < start_time + max_seconds)
    time_left = max_seconds - (now - start_time)

    x, fbest, num_fevals, termination_reason = generating_set_search(f, n; 
      max_seconds = time_left, delta_tol = delta_tol,
      step_size = step_size, x = false, known_fmin = known_fmin, ftol = ftol,
      direction_gen = direction_gen)

    num_runs += 1
    step_size = step_size * 2

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

@time generating_set_search(sphere, 2; known_fmin = 0.0)
@time generating_set_search(rosenbrock, 16; known_fmin = 0.0)
@time restart_gss(rosenbrock, 32; known_fmin = 0.0)

@time restart_gss(rosenbrock, 16; known_fmin = 0.0)
@time restart_gss(xtransform(16, rosenbrock), 16; known_fmin = 0.0)

@time restart_gss(cigar, 32; known_fmin = 0.0)

function repeated_runs(searchf, n = 16, num_runs = 10)
  fevals = zeros(num_runs)
  fbests = zeros(num_runs)
  xs = Any[]
  for(i in 1:num_runs)
    x, fbests[i], fevals[i], termination_reason = searchf(n)
    push!(xs, x)
  end
  return fbests, fevals
end

n = 8
@time fbs, fevals = repeated_runs(
  ((n) -> restart_gss(xtransform(n, rosenbrock), n; known_fmin = 0.0)),
  n)