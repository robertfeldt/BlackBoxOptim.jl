function rosenbrock(x)
  n = length(x)
  return( sum( 100*( x[2:n] - x[1:(n-1)].^2 ).^2 + ( x[1:(n-1)] - 1 ).^2 ) )
end

function rastrigin(x)
  D = length(x)
  10 * D + sum( x.^2 ) - 10 * sum( cos( 2 * π * x ) )
end

function ackley(x)
  D = length(x)
  try
    20 - 20.*exp(-0.2.*sqrt(sum(x.^2)/D)) - exp(sum(cos(2 * π * x))/D) + e
  catch e
    # Sometimes we have gotten a DomainError from the cos function so we protect this call
    println(e)
    println("For input x = ", x)
    # Return a very large fitness value to indicate that this is NOT the solution we want.
    # TODO: Fix better way to handle this!
    9.99e100
  end
end

# Ensure it has been compiled
x, f, es = adaptive_coordinate_descent(rosenbrock, 2, -5.0, 5.0; stopfitness = 1e-10, howOftenUpdateRotation = 1)

function runexp(P, hs, reps = 3)
  nhs = length(hs)

  fs = zeros(nhs, reps)
  fevals = zeros(nhs, reps)
  times = zeros(nhs, reps)

  for hi in 1:nhs
    h = hs[hi]
    for i in 1:reps
      tic()
      xmean, fs[hi, i], fevals[hi, i] = adaptive_coordinate_descent(rosenbrock, P, -5.0, 5.0;
        stopfitness = 1e-10, howOftenUpdateRotation = h)
      times[hi, i] = toc()
    end
  end

  return (mean(fs, 2), mean(fevals, 2), mean(times, 2))
end

# runexp(5, [1, 2, 4, 8, 16], 100) # For P = 5, h = 4-8 seem best speedwise
# runexp(10, [1, 2, 4, 8, 16], 25) # For P = 10, h = 2-8 seem best
# runexp(100, [1, 2, 4, 8, 16], 3) # For P = 100, h = 2-4 seem best

# ACD seems to have problems with schwefel2_21 in 64-dim. rastrigin and ackley is also hard.

func = ackley
P = 16
x, f, es = adaptive_coordinate_descent(func, P, -10.0, 10.0; 
  stopfitness = 1e-10, howOftenUpdateRotation = 2,
  nonImprovementBudget = 1)
