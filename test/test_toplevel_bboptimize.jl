@testset "Top-level interface" begin
    rosenbrock(x) = sum(i -> 100*abs2(x[i+1] - x[i]^2) + abs2(x[i] - 1), Base.OneTo(length(x)-1))

    # Ensure things have compiled before we start testing:
    bboptimize(rosenbrock; SearchRange = (-5.0, 5.0), NumDimensions = 2,
        MaxSteps = 10, TraceMode = :silent)


    @testset "run a simple optimization" begin
        @testset "using bboptimize() with mostly defaults" begin
            res = bboptimize(rosenbrock; SearchRange = (-5.0, 5.0), NumDimensions = 2,
                MaxSteps = 2000, TraceMode = :silent)
            @test best_fitness(res) < 0.1
            xbest = best_candidate(res)
            @test typeof(xbest) == Vector{Float64}
            @test length(xbest) == 2

            # We also mimic some of the Optim.jl api (although it is a bit strange...)
            @test f_minimum(res) < 0.1
            @test minimum(res) == xbest
            @test iteration_converged(res)
        end

        @testset "using bbsetup()/bboptimize() with mostly defaults" begin
            opt = bbsetup(rosenbrock; SearchRange = (-5.0, 5.0), NumDimensions = 2,
                MaxSteps = 2000, TraceMode = :silent)
            @test numruns(opt) == 0
            @test isa(problem(opt), BlackBoxOptim.FunctionBasedProblem)
            res = bboptimize(opt)
            @test numruns(opt) == 1
            @test isa(lastrun(opt), BlackBoxOptim.OptRunController)

            @test best_fitness(res) < 0.1
            xbest = best_candidate(res)
        end

        @testset "using non-population optimizer" begin
            res = bboptimize(rosenbrock; Method=:generating_set_search,
                                            SearchRange = (-5.0, 5.0), NumDimensions = 2,
                                            MaxSteps = 5000, TraceMode = :silent)
            @test best_fitness(res) < 1.0
            xbest = best_candidate(res)
            @test typeof(xbest) == Vector{Float64}
        end

        @testset "using population optimizer" begin
            res = bboptimize(rosenbrock; Method=:adaptive_de_rand_1_bin,
                                            SearchRange = (-5.0, 5.0), NumDimensions = 2,
                                            MaxSteps = 5000, TraceMode = :silent)
            @test best_fitness(res) < 0.1
            xbest = best_candidate(res)
            @test typeof(xbest) == Vector{Float64}
            xpop = population(res)
            @test isa(xpop, BlackBoxOptim.Population)
            @test popsize(xpop) > 0
            @test numdims(xpop) == 2

            # Access a few individual solution vectors in the population...
            @test isa(xpop[1], Array{Float64, 1}) == true # first solution vector
            @test isa(xpop[popsize(xpop)], Array{Float64, 1}) == true # last solution vector
            rand_solution_idx = rand(1:popsize(xpop))
            @test isa(xpop[rand_solution_idx], Array{Float64, 1}) == true # random solution vector

            # and to access their fitness values:
            @test isa(fitness(xpop, 1), Float64)
            @test isa(fitness(xpop, popsize(xpop)), Float64)
            @test isa(fitness(xpop, rand_solution_idx), Float64)

            # Ensure the lowest fitness value is the one returned by best_fitness
            min_fitness_value = minimum(map(i -> fitness(xpop, i), 1:popsize(xpop)))
            @test (min_fitness_value == best_fitness(res))
        end

        @testset "using population optimizer and parallel evaluator" begin
            opt = bbsetup(rosenbrock; Method=:adaptive_de_rand_1_bin,
                                        SearchRange = (-5.0, 5.0), NumDimensions = 2,
                                        MaxSteps = 2000, TraceMode = :silent, Workers=workers())
            res = bboptimize(opt)
            @test isa(BlackBoxOptim.evaluator(lastrun(opt)), BlackBoxOptim.ParallelEvaluator)
        end

    end

    @testset "continue running an optimization after it finished" begin
        optctrl = bbsetup(rosenbrock; SearchRange = (-5.0, 5.0), NumDimensions = 100,
            MaxTime = 0.5, TraceMode = :silent)

        res1 = bboptimize(optctrl)
        @test numruns(optctrl) == 1

        res2 = bboptimize(optctrl; MaxTime = 1.0)
        @test numruns(optctrl) == 2

        @test best_fitness(res2) <= best_fitness(res1)

        # parameters should be the same except for MaxTime
        for p in keys(flatten(parameters(res1)))
            if p != :MaxTime
                @test parameters(res1)[p] == parameters(res2)[p]
            end
        end
        @test parameters(res1)[:MaxTime] == 0.5
        @test parameters(res2)[:MaxTime] == 1.0
    end

    @testset "fitness decrease monotonically if optimizing with same optctrl repeatedly" begin
        # There was a bug in 0.3.0 that did not give monotonically decreasing
        # fitness with repeated runs. The bug only happened for very short MaxSteps
        # lengths, such as 10. This test ensures it does not resurface.
        optctrl = bbsetup(rosenbrock; SearchRange = (-5.0, 5.0), NumDimensions = 100,
            MaxSteps = 10, TraceMode = :silent)
        optresults = [bboptimize(optctrl) for _ in 1:10]
        bestfitnesses = map(best_fitness, optresults)
        for i in 2:length(optresults)
            @test bestfitnesses[i-1] >= bestfitnesses[i]
        end
    end

    @testset "return results after interruption" begin
        i = 0
        function rosenbrock_throwing(x)
            i += 1
            if i < 50
                return sum(i -> 100*(x[i+1] - x[i])^2 + (x[i] - 1)^2, 1:(length(x)-1))
            else
                throw(InterruptException())
            end
        end
        @testset ":RecoverResults on" begin
            i = 0
            optctrl = bbsetup(rosenbrock_throwing; SearchRange = (-5.0, 5.0), NumDimensions = 100,
                              MaxSteps=100, TraceMode=:silent, RecoverResults=true)
            res = bboptimize(optctrl)
            @test BlackBoxOptim.isinterrupted(res)
        end

        @testset ":RecoverResults off" begin
            i = 0 # reset the counter, otherwise it will throw in the setup
            optctrl = bbsetup(rosenbrock_throwing; SearchRange = (-5.0, 5.0), NumDimensions = 100,
                              MaxSteps=100, TraceMode=:silent, RecoverResults=false)
            @test_throws InterruptException bboptimize(optctrl)
        end
    end

    @testset "continue running an optimization after serializing to disc" begin
        using Serialization: serialize, deserialize

        optctrl = bbsetup(rosenbrock; SearchRange = (-5.0, 5.0), NumDimensions = 100,
                          MaxTime = 0.5, TraceMode = :silent)
        res1 = bboptimize(optctrl)

        local tempfilename = "./temp" * string(rand(1:10^8)) * ".tmp"

        try # To ensure we delete the temp file afterwards...
            open(tempfilename, "w") do fh
                serialize(fh, optctrl)
            end

            # Try to make sure its not in mem:
            opctrl = nothing; GC.gc()

            local ocloaded
            open(tempfilename, "r") do fh
                ocloaded = deserialize(fh)
            end

            @test numruns(ocloaded) == 1
            res2 = bboptimize(ocloaded; MaxTime = 1.0)
            @test numruns(ocloaded) == 2

            @test best_fitness(res2) <= best_fitness(res1)
        finally
            if isfile(tempfilename)
                rm(tempfilename)
            end
        end
    end

    @testset "TargetFitness option works" begin
        # FIXME use the same (fixed?) random seed to guarantee reproducibility
        result1 = bboptimize(rosenbrock, SearchRange = (-5.0, 5.0), NumDimensions = 5,
                             Method = :de_rand_1_bin, FitnessTolerance = 1e-5,
                             MaxSteps = 1000000, TraceMode = :silent,
                             TargetFitness = 0.0)
        result2 = bboptimize(rosenbrock, SearchRange = (-5.0, 5.0), NumDimensions = 5,
                             Method = :de_rand_1_bin, FitnessTolerance = 1e-5,
                             MaxSteps = 1000000, TraceMode = :silent)
        @test best_fitness(result1) < 1e-5
        @test result1.iterations < result2.iterations
    end

    @testset "callback" begin
        global callbacktimes = Float64[]
        function callbackfn(optctrl)
            global callbacktimes
            push!(callbacktimes, time())
        end
        interval = 0.005
        NumCalls = 30
        res = bboptimize(rosenbrock; SearchRange = (-5.0, 5.0), NumDimensions = 20,
            TraceMode = :silent, 
            MaxTime = 1.0 + interval * NumCalls, # Some extra time for startup/compilation etc
            CallbackInterval = interval, CallbackFunction = callbackfn)
        @test length(callbacktimes) >= NumCalls
        for i in 1:(length(callbacktimes) - 1)
            d = (callbacktimes[i+1] - callbacktimes[i])
            @test d >= interval
        end

        # If we call with an interval which is negative there are no callbacks
        prenumcalls = length(callbacktimes)
        res = bboptimize(rosenbrock; SearchRange = (-5.0, 5.0), NumDimensions = 20,
            TraceMode = :silent, 
            MaxTime = 1.0 + interval * NumCalls, # Some extra time for startup/compilation etc
            CallbackInterval = -2.0, CallbackFunction = callbackfn)
        @test length(callbacktimes) == prenumcalls
    end
end
