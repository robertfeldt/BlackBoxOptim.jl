using BlackBoxOptim, HTTP, Sockets

# Create the fitness plot object and serve it
const vlplot = BlackBoxOptim.VegaLiteMetricOverTimePlot(verbose=false)
HTTP.serve!(vlplot)

# Func to optimize.
function rosenbrock(x)
    sum(i -> 100*abs2(x[i+1] - x[i]^2) + abs2(x[i] - 1), Base.OneTo(length(x)-1))
end

# Now optimize for 2 minutes.
# Go to http://127.0.0.1:8081 to view fitness progress!
res = bboptimize(rosenbrock;
        SearchRange=(-10.0,10.0), NumDimensions = 500,
        PopulationSize=100, MaxTime=2*60.0,
        CallbackFunction = Base.Fix1(BlackBoxOptim.fitness_plot_callback, vlplot),
        CallbackInterval = 2.0);
println("Best fitness = ", best_fitness(res))
