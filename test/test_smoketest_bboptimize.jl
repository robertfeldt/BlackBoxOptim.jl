@testset "bboptimize() single-objective methods smoketest" begin
    rosenbrock2d(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2

    for m in keys(BlackBoxOptim.SingleObjectiveMethods)
        @testset "$(m)" begin
            ctrl = bbsetup(rosenbrock2d; Method = m,
                SearchRange = [(-5.0, 5.0), (-2.0, 2.0)], TraceMode = :silent)
            # run first iteration before the main run to exclude compilation from timing
            bboptimize(ctrl, MaxSteps = 1)
            res = bboptimize(ctrl, MaxTime = 0.3)
            @test length(best_candidate(res)) == 2
            f = best_fitness(res)
            @test typeof(f) == Float64
            @test f < 100.0 # this can't be very tight since we give very little time for optimization...
        end
    end
end

@testset "bboptimize() multi-objective methods smoketest" begin
    schaffer1(x) = (sumabs2(x), sumabs2(x .- 2.0))

    for m in keys(BlackBoxOptim.MultiObjectiveMethods)
        @testset "$(m)" begin
            ctrl = bbsetup(schaffer1; Method = m, 
                SearchRange = [(-10.0, 10.0), (-10.0, 10.0)], TraceMode = :silent,
                FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true), ϵ=0.01)
            # run first iteration before the main run to exclude compilation from timing
            bboptimize(ctrl, MaxSteps = 1)
            res = bboptimize(ctrl, MaxTime = 3.0)
            @test length(best_candidate(res)) == 2
            f = best_fitness(res)
            @test typeof(f) == NTuple{2,Float64}
            @test isapprox(f[1], 2.0; atol=8E-2)
            @test isapprox(f[2], 2.0; atol=8E-2)
        end
    end
end
