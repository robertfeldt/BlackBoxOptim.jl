using BlackBoxOptim

function rosenbrock(x)
  return( sum( 100*( x[2:end] - x[1:end-1].^2 ).^2 + ( x[1:end-1] - 1 ).^2 ) )
end

# First run without any parallel procs used in eval
opt1 = bbsetup(rosenbrock; SearchRange = (-5.0, 5.0),
  NumDimensions = 100, MaxSteps = 500000)
res1 = bboptimize(opt1)

# Now add 2 procs that can exec in parallel (obviously it depends on your CPU
# what you actually gain from this though)
addprocs(2)
# You need to make sure BlackBoxOptim is loaded an all procs:
@everywhere using BlackBoxOptim
# Now we can setup a new opt and run the optimization:
opt2 = bbsetup(rosenbrock; SearchRange = (-5.0, 5.0),
  NumDimensions = 100, MaxSteps = 500000, Workers = workers())
res2 = bboptimize(opt2)
