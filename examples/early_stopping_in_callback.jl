using BlackBoxOptim

rosenbrock(x) = sum( 100*( x[2:end] .- x[1:end-1].^2 ).^2 .+ ( x[1:end-1] .- 1 ).^2 )

const MyFitnessGoal = 30.0

function myfitnessgoalachieved(oc)
    best_fitness(oc) < MyFitnessGoal
end

function cbearlystopping(oc)
    if myfitnessgoalachieved(oc)
        BlackBoxOptim.shutdown!(oc)
    end
end

# We give max 1000 seconds to optimize but expect this to stop much earlier
@time res = bboptimize(rosenbrock; SearchRange = (-100.0, 100.0), NumDimensions = 100, 
        CallbackFunction = cbearlystopping, CallbackInterval = 0.0, MaxTime = 1000.0);

@assert best_fitness(res) < MyFitnessGoal
@assert BlackBoxOptim.stop_reason(res) == "Run explicitly stopped via shutdown method"