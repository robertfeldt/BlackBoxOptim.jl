using Distributed

# Now add 2 procs that can exec in parallel (obviously it depends on your CPU
# what you actually gain from this though)
addprocs(2)

# Ensure BlackBoxOptim loaded on all workers
@everywhere using BlackBoxOptim

# define the function to optimize on all workers. Parallel eval only gives a gain
# if function to optimize is slow. For this example we introduce a fake sleep
# to make it slow since the function is actually very quick to eval...
@everywhere function slow_rosenbrock(x)
  sleep(0.001) # Fake a slower func to be optimized...
  return BlackBoxOptim.rosenbrock(x)
end

# First run without any parallel procs used in eval
opt1 = bbsetup(slow_rosenbrock; Method=:xnes, SearchRange = (-5.0, 5.0),
               NumDimensions = 50, MaxFuncEvals = 5000)
el1 = @elapsed res1 = bboptimize(opt1)
t1 = round(el1, digits=3)

# When Workers= option is given, BlackBoxOptim enables parallel
# evaluation of fitness using the specified worker processes
opt2 = bbsetup(slow_rosenbrock; Method=:xnes, SearchRange = (-5.0, 5.0),
               NumDimensions = 50, MaxFuncEvals = 5000, Workers = workers())
el2 = @elapsed res2 = bboptimize(opt2)
t2 = round(el2, digits=3)

println("Time: serial = $(t1)s, parallel = $(t2)s")
if t2 < t1
  println("Speedup is $(round(t1/t2, digits=1))x")
else
  println("Slowdown is $(round(t2/t1, digits=1))x")
end
