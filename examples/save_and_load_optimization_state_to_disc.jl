using BlackBoxOptim
using Serialization

# Lets start an optimization of 2D rosenbrock
# as in the README, then save its state to disc,
# read it back in and continue optimization.
rosenbrock2d(x) = abs2(1.0 - x[1]) + 100.0 * abs2(x[2] - x[1]^2)

# If you want to restart optimization you need to get a handle
# to an optimization controller. Use the function bbsetup
# instead of bboptimize, but with the same parameters.
# We just run it for 100 steps though:
optctrl = bbsetup(rosenbrock2d; SearchRange = (-5.0, 5.0), NumDimensions = 2,
    MaxSteps = 100, TraceMode = :silent);

# Now we can run it:
res100 = bboptimize(optctrl)

# Print the best and a randomly selected candidate solution so we can ensure
# they are the same later.
best100  = best_candidate(res100)
idx = rand(1:popsize(optctrl.optimizer.population))
acand100 = optctrl.optimizer.population[idx]
println("Best candidate: ", best100)
println("Candidate num $(idx): ", acand100)

# Now serialize to a temp file:
tempfilename = "./temp" * string(rand(1:Int(1e8))) * ".tmp"
fh = open(tempfilename, "w")
serialize(fh, (optctrl, res100))
close(fh)

# Read back in from file:
fh = open(tempfilename, "r")
optctrlb, res100b = deserialize(fh);
close(fh)

# Print the same candidates:
best100b  = best_candidate(res100b)
acand100b = optctrlb.optimizer.population[idx]
println("Best candidate after load: ", best100)
println("Candidate num $(idx) after load: ", acand100)

# Clean up the temp file:
rm(tempfilename)

# Now restart the optimization (but run more steps):
res1100 = bboptimize(optctrlb; MaxSteps = 1000)

# And print the fitness progress:
println("Fitness progress: ", (best_fitness(res100), best_fitness(res1100)))

# NOTE! This method is not guaranteed to give stable results when you change
# version of julia between saving and loading the data.
