using BlackBoxOptim

# Lets start and then restart an optimization of 2D rosenbrock
# as in the README:
rosenbrock2d(x) = abs2(1.0 - x[1]) + 100.0 * abs2(x[2] - x[1]^2)

# If you want to restart optimization you need to get a handle
# to the optimizer, problem and params. Use the function
# setup_bboptimize instead of bboptimize, but with the same
# parameters. We just run it for 10 steps though:
opt = bbsetup(rosenbrock2d; SearchRange = (-5.0, 5.0), NumDimensions = 2,
    MaxSteps = 10, TraceMode = :silent);

# Now we can run it (this will do 10 steps):
res10 = bboptimize(opt)

# And lets run it again (another 10 steps):
res20 = bboptimize(opt)

# Lets run 980 steps this time to try to get more progress:
res1000 = bboptimize(opt; MaxSteps = 980)

# And finish it off with another 9000 steps in batches of 1000:
results = [bboptimize(opt; MaxSteps = 1000) for _ in 1:9]
res10k = results[end]

# Now compare fitness progress:
println("Fitness progress: ", map(best_fitness, [res10, res20, res1000, res10k]))

# And note that the elapsed time was much larger for first call since things
# were compiled and that the elapsed time for last two entries should be about the same
# (since they were with 980 and 1000 steps, respectively):
println("Elapsed time: ", map(BlackBoxOptim.elapsed_time, [res10, res20, res1000, res10k]))

# And results should approach the Rosenbrock optimum at [1.0, 1.0]:
println("Best solution progress: ", map(best_candidate, [res10, res20, res1000, res10k]))
