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
function rosenbrock(x)
    sum(i -> 100*abs2(x[i+1] - x[i]^2) + abs2(x[i] - 1), Base.OneTo(length(x)-1))
end

# Now optimize for 2 minutes. 
# Go to http://127.0.0.1:8081 to view fitness progress!
res = bboptimize(rosenbrock;
        SearchRange=(-10.0,10.0), NumDimensions = 500,
        PopulationSize=100, MaxTime=2*60.0,
        CallbackFunction = callback, CallbackInterval = 2.0);
println("Best fitness = ", best_fitness(res))