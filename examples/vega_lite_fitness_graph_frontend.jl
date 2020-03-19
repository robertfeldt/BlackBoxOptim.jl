using BlackBoxOptim

# Create the fitness graph object and serve it
const Vlg = BlackBoxOptim.VegaLiteGraphFitnessGraph(false)
BlackBoxOptim.serve(Vlg)

# We will use a callback function to add fitness data to the graph.
# Since no "Time" is given the graph object will set time 0.0 for the first
# time we add data and then will count relative to that.
function callback(oc)
    BlackBoxOptim.add_data!(Vlg, Dict("Fitness" => best_fitness(oc)))
end

# Func to optimize.
function rastrigin(x)
    D = length(x)
    10 * D + sum(abs2, x) - 10 * sum(xx -> cos(2π * xx), x)
end

# Now optimize for 2 minutes. 
# Go to http://127.0.0.1:8081 to view fitness progress!
res = bboptimize(rastrigin;
        SearchRange=(-5.12,5.12), NumDimensions = 50,
        PopulationSize=100, MaxTime=2*60.0,
        CallbackFunction = callback, CallbackInterval = 2.0);
println("Best fitness = ", best_fitness(res))