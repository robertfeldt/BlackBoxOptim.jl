using BlackBoxOptim
using BlackBoxOptim: num_func_evals # We will update progress based on num func evals so far
using ProgressMeter

# Case 1: Optimization for a given number of function evaluations
MaxFuncEvals = 5000000

# Minimum update interval 0.5 seconds
const Prog = Progress(MaxFuncEvals, 0.5, "Optimizing...") 

function callback_progress_stepper(optcontroller)
    global Prog
    ProgressMeter.update!(Prog, num_func_evals(optcontroller))
end

# Let's try this on 100-dimensional rosenbrock
function rosenbrock(x)
    return( sum( 100*( x[2:end] .- x[1:end-1].^2 ).^2 .+ ( x[1:end-1] .- 1 ).^2 ) )
end

println("Starting optimization!")
res = bboptimize(rosenbrock;
        SearchRange=(-5.0,5.0), NumDimensions = 100,
        PopulationSize=10, MaxFuncEvals=MaxFuncEvals,
        CallbackFunction = callback_progress_stepper, 
        CallbackInterval = 0.0,
        TraceMode = :silent # So that output from BBO doesn't interfere with progress bar...
);

# We can adapt this to other use cases (such as max time etc) by keeping track of start time
# and then updating based on delta time etc.