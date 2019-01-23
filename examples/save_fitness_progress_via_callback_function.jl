using BlackBoxOptim

# We also need the (not yet export) num_func_evals
using BlackBoxOptim: num_func_evals

function rosenbrock2d(x)
    return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
end

fitness_progress_history = Array{Tuple{Int, Float64},1}()
callback = oc -> push!(fitness_progress_history, (num_func_evals(oc), best_fitness(oc)))

# A callback interval of 0.0 ensures we callback every step of the optimization.
res = bboptimize(rosenbrock2d; 
        SearchRange=(-5.0,5.0), NumDimensions = 2, 
        PopulationSize=10, MaxFuncEvals=200, 
        CallbackFunction = callback, CallbackInterval = 0.0);

# Note that fitness_progress_history is NOT guaranteed to have all values for the 
# number of function evaluations since it depends on the optimization algorithm how 
# many function evaluations they do in each step. The callback happen only between 
# steps. The number of func evals might even be stochastic or change per step.
println(fitness_progress_history)