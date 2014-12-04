using BlackBoxOptim

# Lets start an optimization of 2D rosenbrock
# as in the README, then save its state to disc, 
# read it back in and continue optimization.
function rosenbrock2d(x)
  return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
end

# If you want to restart optimization you need to get a handle
# to the optimizer, problem and params. Use the function
# setup_bboptimize instead of bboptimize, but with the same
# parameters. We just run it for 10 steps though:
optimizer, problem, params = BlackBoxOptim.setup_bboptimize(rosenbrock2d; 
  search_range = (-5.0, 5.0), dimensions = 2, 
  parameters = {:MaxSteps => 100, :ShowTrace => false});

# Now we can run it:
b100, f100, tr100, time100, params, ne100 = BlackBoxOptim.run_optimizer_on_problem(
  optimizer, problem; parameters = params);

# Now serialize to a temp file:
tempfilename = "./temp" * string(rand(1:int(1e8))) * ".tmp"
fh = open(tempfilename, "w")
serialize(fh, (optimizer, problem, params))
close(fh)

# Read back in from file:
fh = open(tempfilename, "r")
opt2, prob2, params2 = deserialize(fh)
close(fh)

# Clean up the temp file:
rm(tempfilename)

# Now restart the optimization:
b200, f200, tr200, time200, params, ne200 = BlackBoxOptim.run_optimizer_on_problem(
  opt2, prob2; parameters = params2);

# And print the fitness progress:
println("Fitness progress: ", (f100, f200))

# NOTE! This method is not guaranteed to give stable results when you change
# version of julia between saving and loading the data.