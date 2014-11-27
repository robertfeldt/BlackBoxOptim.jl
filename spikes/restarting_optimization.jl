using BlackBoxOptim

# Lets start and then restart an optimization of 2D rosenbrock
# as in the README:
function rosenbrock2d(x)
  return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
end

# If you want to restart optimization you need to get a handle
# to the optimizer, problem and params. Use the function
# setup_bboptimize instead of bboptimize, but with the same
# parameters. We just run it for 10 steps though:
optimizer, problem, params = BlackBoxOptim.setup_bboptimize(rosenbrock2d; 
  search_range = (-5.0, 5.0), dimensions = 2, 
  parameters = {:MaxSteps => 10, :ShowTrace => false});

# Now we can run it:
best10, fitness10, termination_reason10, elapsed_time10, params, num_evals10 = BlackBoxOptim.run_optimizer_on_problem(optimizer, problem; parameters = params);

# And lets run it again:
best20, fitness20, termination_reason20, elapsed_time20, params, num_evals20 = BlackBoxOptim.run_optimizer_on_problem(optimizer, problem; parameters = params);

# Lets run 980 steps this time to try to get more progress:
params[:MaxSteps] = 980
best1000, fitness1000, termination_reason1000, elapsed_time1000, params, num_evals1000 = BlackBoxOptim.run_optimizer_on_problem(optimizer, problem; parameters = params);

# And finish it off with another 9000 steps in batches of 100 (and without tracing):
params[:MaxSteps] = 1000
fitness10000 = best10000 = elapsed_time1000b = 1 # Just so saved outside of loop body...
for i in 1:9
  best10000, fitness10000, termination_reason10000, elapsed_time1000b, params, num_evals10000 = BlackBoxOptim.run_optimizer_on_problem(optimizer, problem; parameters = params);
end

# Now compare fitness progress:
println("Fitness progress: ", (fitness10, fitness20, fitness1000, fitness10000))

# And note that the elapsed time was much larger for first call since things
# were compiled and that the elapsed time for last two entries should be about the same
# (since they were with 980 and 1000 steps, respectively):
println("Elapsed time: ", (elapsed_time10, elapsed_time20, elapsed_time1000, elapsed_time1000b))

# And results should approach the Rosenbrock optimum at [1.0, 1.0]:
println("Best solution progress: ", (best10, best20, best1000, best10000))