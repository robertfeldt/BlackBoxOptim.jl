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

# Average number of func evaluations to reach a target value divided by success rate.
function bbob_sp1(func, targetValue, P, nonImprovementBudget = 1, reps = 10, MaxEval = 1e4*P)
  nfevals = zeros(Int, reps)
  fitnesses = zeros(reps)
  times = zeros(reps)
  for r in 1:reps
    tic()
    x, fitnesses[r], nfevals[r] = adaptive_coordinate_descent(func, P, -10.0, 10.0; 
      MaxEval = MaxEval,
      stopfitness = targetValue, howOftenUpdateRotation = 1,
      nonImprovementBudget = 1)
    times[r] = toq()
  end
  successes = fitnesses .<= targetValue
  numsuccesses = sum(successes)
  successrate = numsuccesses / reps
  mean_fevals_successes = mean(nfevals[successes])
  mean_times_successes = mean(times[successes])
  sp1 = mean_fevals_successes / successrate
  std_fevals = round(std(nfevals[successes]), digits=2)
  std_times = round(std(times[successes]), digits=2)

  println("\n\nNum fevals: $(round(mean_fevals_successes, digits=2)) +/- $(std_fevals) ($(minimum(nfevals[successes]))-$(maximum(nfevals[successes])))")
  println("Times     : $(round(mean_times_successes, digits=2)) +/- $(std_times) ($(minimum(times[successes]))-$(maximum(times[successes])))")
  println("SP1 = $(round(sp1, digits=2))")
  println("log10(fevals/dim) = $(round(log10(mean_fevals_successes/P), digits=2))")

  return (sp1, successrate)
end

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
      times[hi, i] = toq()
    end
  end

  return (mean(fs, 2), mean(fevals, 2), mean(times, 2))
end

# Ensure it has been compiled
x, f, es = adaptive_coordinate_descent(rosenbrock, 2, -5.0, 5.0; stopfitness = 1e-10, howOftenUpdateRotation = 1)
x, f, es = adaptive_coordinate_descent(rastrigin, 10, -100.0, 100.0; stopfitness = 1e-8, howOftenUpdateRotation = 1)

# runexp(5, [1, 2, 4, 8, 16], 100) # For P = 5, h = 4-8 seem best speedwise
# runexp(10, [1, 2, 4, 8, 16], 25) # For P = 10, h = 2-8 seem best
# runexp(100, [1, 2, 4, 8, 16], 3) # For P = 100, h = 2-4 seem best

# ACD seems to have problems with schwefel2_21 in 64-dim. rastrigin and ackley is also hard.

fitnessfct = rosenbrock
P = 2
Xmin = -5.0
Xmax = 5.0
MaxEval = 1e4 * P
stopfitness = 1e-8
howOftenUpdateRotation = 1
numfevals = 0
kSuccess = 2.0
kUnsuccess = 0.5
alwaysSampleTwoPoints = false
nonImprovementBudget = 1

#x, f, es = adaptive_coordinate_descent(fitnessfct, P, -10.0, 10.0; 
#  stopfitness = 1e-8, howOftenUpdateRotation = P,
#  nonImprovementBudget = 1)

# Cmu = 0, MaxEval = 1e4*P.
bbob_sp1(rosenbrock, 1e-8, 2, 1, 30)  # 422 +/- 184
bbob_sp1(rosenbrock, 1e-8, 10, 1, 30) # 5081 +/- 2226
bbob_sp1(rosenbrock, 1e-8, 30, 1, 30) # 49568 +/- 21117

P = 100
bbob_sp1(rosenbrock, 1e-8, P, 1, 3, 1e5*P) # 978986.0 +/- 397922.39 (578444-1374236)
bbob_sp1(rosenbrock, 1e-8, P, iceil(sqrt(P)), 3, 1e5*P) # 1.68466367e6 +/- 566216.47 (1322155-2337131)
bbob_sp1(rosenbrock, 1e-8, P, 4, 3, 1e5*P) # 978986.0 +/- 397922.39 (578444-1374236)

P = 256
bbob_sp1(rosenbrock, 1e-8, P, iceil(sqrt(P)), 3, 1e5*P) # 

bbob_sp1(ackley, 1e-4, 10, 1, 30) # 
bbob_sp1(ackley, 1e-4, 30, 1, 10) # 
