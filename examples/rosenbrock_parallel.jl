# Now add 2 procs that can exec in parallel (obviously it depends on your CPU
# what you actually gain from this though)
addprocs(2)

# note that BlackBoxOptim gets automatically used on all workers
using BlackBoxOptim

# define the function to optimize on all workers. Parallel eval only gives a gain
# if function to optimize is slow. For this example we introduce a fake sleep
# to make it slow since the function is actually very quick to eval...
@everywhere function rosenbrock(x)
  sleep(0.001) # Fake a slower func to be optimized...
  return( sum( 100*( x[2:end] - x[1:end-1].^2 ).^2 + ( x[1:end-1] - 1 ).^2 ) )
end

# First run without any parallel procs used in eval
opt1 = bbsetup(rosenbrock; Method=:xnes, SearchRange = (-5.0, 5.0),
               NumDimensions = 50, MaxFuncEvals = 5000)
tic()
res1 = bboptimize(opt1)
t1 = round(toq(), 3)

# When Workers= option is given, BlackBoxOptim enables parallel
# evaluation of fitness using the specified worker processes
opt2 = bbsetup(rosenbrock; Method=:xnes, SearchRange = (-5.0, 5.0),
               NumDimensions = 50, MaxFuncEvals = 5000, Workers = workers())
tic()
res2 = bboptimize(opt2)
t2 = round(toq(), 3)

println("Time: serial = $(t1)s, parallel = $(t2)s")
if t2 < t1
  println("Speedup is $(round(t1/t2, 1))x")
else
  println("Slowdown is $(round(t2/t1, 1))x")
end
