# Now add 2 procs that can exec in parallel (obviously it depends on your CPU
# what you actually gain from this though)
addprocs(2)

# note that BlackBoxOptim gets automatically used on all workers
using BlackBoxOptim

# define the function to optimize on all workers
@everywhere function rosenbrock(x)
  return( sum( 100*( x[2:end] - x[1:end-1].^2 ).^2 + ( x[1:end-1] - 1 ).^2 ) )
end

# First run without any parallel procs used in eval
opt1 = bbsetup(rosenbrock; Method=:xnes, SearchRange = (-5.0, 5.0),
               NumDimensions = 100, MaxSteps = 5000)
res1 = bboptimize(opt1)

# When Workers= option is given, BlackBoxOptim enables parallel
# evaluation of fitness using the specified worker processes
opt2 = bbsetup(rosenbrock; Method=:xnes, SearchRange = (-5.0, 5.0),
               NumDimensions = 100, MaxSteps = 5000, Workers = workers())
res2 = bboptimize(opt2)
