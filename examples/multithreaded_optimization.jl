using Pkg
Pkg.activate(".")

if Threads.nthreads() < 2
    error("Multiple threads NOT found! Makes no sense to test this without them... Start as, for example: JULIA_NUM_THREADS=4 julia multithreaded_optimization.jl")
    exit(-1)
end

using BlackBoxOptim

# Functions to optimize. Should be thread-safe.
function rosenbrock(x)
    sleep(0.05) # So that there is some benefit in the thread switching...
    sum(i -> 100*abs2(x[i+1] - x[i]^2) + abs2(x[i] - 1), Base.OneTo(length(x)-1))
end

function rastrigin(x)
    D = length(x)
    10 * D + sum(abs2, x) - 10 * sum(xx -> cos(2π * xx), x)
end

# For multi-objective opt:
optfun(x) = (rosenbrock(x), rastrigin(x))

# Let's ensure they are really thread-safe:
for (fn, N) in [(rosenbrock, 100), (rastrigin, 1000000)]
    points = rand(N, 2)
    @time single_threaded = map(i -> fn(points[i, :]), 1:N)
    multi_threaded = zeros(Float64, N)
    @time Threads.@threads for i in 1:N
        multi_threaded[i] = fn(points[i, :])
    end
    @assert sum((single_threaded .- multi_threaded).^2) == 0
end

# Short runs to ensure things have compiled:
bboptimize(rosenbrock; SearchRange = (-5.0, 5.0), NumDimensions = 2,
    MaxTime = 0.1, TraceMode = :silent, Method = :dxnes, lambda = 10);
bboptimize(rosenbrock; SearchRange = (-5.0, 5.0), NumDimensions = 2,
    MaxTime = 0.1, TraceMode = :silent, NThreads=Threads.nthreads()-1, Method = :dxnes, lambda = 10);
bboptimize(optfun; SearchRange = (-5.0, 5.0), NumDimensions = 2,
    MaxTime = 0.5, TraceMode = :silent, NThreads=Threads.nthreads()-1,
    Method = :borg_moea, ϵ=0.05,
    FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true));

# Same time given to single- and multi-threaded optimization runs:
MaxTime = 30.0

# Now run the single-threaded optimization:
@time res_single = bboptimize(rosenbrock; SearchRange = (-5.0, 5.0), NumDimensions = 100,
    MaxTime = MaxTime, TraceMode = :silent, Method = :dxnes, lambda = 50)

# Now run the multi-threaded optimization:
@time res_multi = bboptimize(rosenbrock; SearchRange = (-5.0, 5.0), NumDimensions = 100,
    MaxTime = MaxTime, TraceMode = :silent, NThreads=Threads.nthreads()-1,
    Method = :dxnes, lambda = 50)

println("Fitness (single-threaded): ", round(best_fitness(res_single), digits=4))
println("Fitness (multi-threaded):  ", round(best_fitness(res_multi), digits=4))
#@assert BlackBoxOptim.f_calls(res_single) < BlackBoxOptim.f_calls(res_multi)

@time res_single = bboptimize(optfun; SearchRange = (-5.0, 5.0), NumDimensions = 100,
    MaxTime = MaxTime, TraceMode = :silent,
    Method = :borg_moea, ϵ=0.05,
    FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true));

@time res_multi = bboptimize(optfun; SearchRange = (-5.0, 5.0), NumDimensions = 100,
    MaxTime = MaxTime, TraceMode = :silent, NThreads=Threads.nthreads()-1,
    Method = :borg_moea, ϵ=0.05,
    FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true));

# The whole idea with multi-threading is that we can evaluate more function calls in same time so:
@assert BlackBoxOptim.f_calls(res_single) < BlackBoxOptim.f_calls(res_multi)

# Now do repeated runs and compare if there is a real difference in fitness
using HypothesisTests
using Statistics

for method in [:adaptive_de_rand_1_bin_radiuslimited, :dxnes]
    fsingle, fmulti = Float64[], Float64[]
    for i in 1:10
        println("Single-threaded run $i")
        res_single = bboptimize(rosenbrock; SearchRange = (-5.0, 5.0), NumDimensions = 100,
            Method = method, MaxTime = MaxTime, TraceMode = :silent, lambda = 50)
        push!(fsingle, best_fitness(res_single))
        println("Fitness (single-threaded): ", round(last(fsingle), digits=2))

        println("Multi-threaded run $i")
        res_multi = bboptimize(rosenbrock; SearchRange = (-5.0, 5.0), NumDimensions = 100,
            Method = method, MaxTime = MaxTime, TraceMode = :silent, lambda = 50,
            NThreads=Threads.nthreads()-1)
        push!(fmulti, best_fitness(res_multi))
        println("Fitness (multi-threaded): ", round(last(fmulti), digits=2))
    end

    pval = pvalue(MannWhitneyUTest(fsingle, fmulti))

    # Seems to give consistently better fitness in single-threaded runs on my 4-core
    # (as well as 2-core (but more expected there)) laptops for adaptive_de_rand_1_bin_radiuslimited.
    # Tried for MaxTime 5, 10, and 60 seconds with same results...
    println("For method: ", method)
    println("Mean Fitness (single-threaded): ", round(mean(fsingle), digits=2))
    println("Mean Fitness (multi-threaded):  ", round(mean(fmulti), digits=2))
    println("Significant difference: ", (pval < 0.01 ? "Yes" : "No"), " ($(round(pval, digits=5)))")
end
