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
  parameters = {:MaxSteps => 10});

# Now we can run it:
best10, fitness10, termination_reason10, elapsed_time10, params, num_evals10 = BlackBoxOptim.run_optimizer_on_problem(optimizer, problem; parameters = params);

# And lets run it again:
best20, fitness20, termination_reason20, elapsed_time20, params, num_evals20 = BlackBoxOptim.run_optimizer_on_problem(optimizer, problem; parameters = params);

# Lets run 880 steps this time to try to get more progress:
params[:MaxSteps] = 980
best1000, fitness1000, termination_reason1000, elapsed_time1000, params, num_evals1000 = BlackBoxOptim.run_optimizer_on_problem(optimizer, problem; parameters = params);

# And finish it off with another 9000 steps:
params[:MaxSteps] = 9000
best10000, fitness10000, termination_reason10000, elapsed_time10000, params, num_evals10000 = BlackBoxOptim.run_optimizer_on_problem(optimizer, problem; parameters = params);

# Now compare fitness progress:
println((fitness10, fitness20, fitness1000, fitness10000))

# And note that the elapsed time was much larger for first call since things
# were compiled:
println((elapsed_time10, elapsed_time20, elapsed_time1000, elapsed_time10000))