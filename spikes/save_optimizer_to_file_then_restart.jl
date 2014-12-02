using BlackBoxOptim
using HDF5, JLD

# Lets start and then restart an optimization of 2D rosenbrock
# as in the README:
function rosenbrock2d(x)
  return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
end

# If you want to save an optimization state you need to get a handle
# to the optimizer, problem and params. Use the function
# setup_bboptimize instead of bboptimize, but with the same
# parameters. We just run it for 10 steps though:
optimizer, problem, params = BlackBoxOptim.setup_bboptimize(rosenbrock2d; 
  search_range = (-5.0, 5.0), dimensions = 2, 
  parameters = {:MaxSteps => 10, :ShowTrace => false});

# Now we can run it:
best10, fitness10, termination_reason10, elapsed_time10, params, num_evals10 = BlackBoxOptim.run_optimizer_on_problem(optimizer, problem; parameters = params);

# And save it with JLD to disc. But can't save pointers! How to handle this???
filename = "./temp" * string(rand(1:int(1e10))) * ".jld"
save(filename, "optimizer", optimizer, "problem", problem, "params", params)

# And load them back:
d = load("/tmp/myfile.jld")
optimizer2 = d["optimizer"]
problem2 = d["problem"]
params2 = d["params"]

# Lets run 9990 steps this time to try to get more progress:
params[:MaxSteps] = 9990
best10000, fitness10000, termination_reason10000, elapsed_time10000, params, num_evals10000 = BlackBoxOptim.run_optimizer_on_problem(optimizer2, problem2; parameters = params2);

