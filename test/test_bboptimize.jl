rosenbrock2d(x) = 100.0*abs2(x[2] - x[1]^2) + abs2(x[1] - 1.0)

function rosenbrock(x)
    @inbounds res = sum(i -> 100.0*abs2(x[i+1] - x[i]^2) + abs2(x[i] - 1.0), 1:length(x)-1)
    return res
end

@testset "bboptimize" begin
    @testset "README examples" begin
        @testset "example #1" begin
            res = bboptimize(rosenbrock2d; SearchRange = (-5.0, 5.0), NumDimensions = 2, TraceMode = :silent)
            @test best_fitness(res) < 0.001
            @test length(best_candidate(res)) == 2
        end

        @testset "example #2" begin
            res = bboptimize(rosenbrock2d; SearchRange = [(-5.0, 5.0), (-2.0, 2.0)], TraceMode = :silent)
            @test best_fitness(res) < 0.001
        end

        @testset "example #3" begin
            res = bboptimize(rosenbrock2d; SearchRange = (-5.0, 5.0), NumDimensions = 2, Method = :de_rand_1_bin, TraceMode = :silent)
            @test best_fitness(res) < 0.001
        end

        @testset "example #4" begin
            res = bboptimize(rosenbrock2d; SearchRange = (-5.0, 5.0), NumDimensions = 2,
                            Method = :random_search, MaxTime = 10.0, TraceMode = :silent)
            @test best_fitness(res) < 0.2
        end

        @testset "example #5" begin
            res = BlackBoxOptim.compare_optimizers(rosenbrock; SearchRange = (-5.0, 5.0), NumDimensions = 3,
                            MaxTime = 2.0, TraceMode = :compact)

            # We at least expect the DE optimizers and DX-NES to come out better than random search
            idx_adaptive_de = findfirst(t -> t[1] == :adaptive_de_rand_1_bin, res)
            idx_random = findfirst(t -> t[1] == :random_search, res)
            @test idx_adaptive_de < idx_random

            idx_dxnes = findfirst(t -> t[1] == :dxnes, res)
            @test idx_dxnes < idx_random
        end
    end

    # run one longer example in case there is problem with the reporting in long runs
    @testset "long runs reporting" begin
        res = bboptimize(rosenbrock2d; SearchRange = (-5.0, 5.0), NumDimensions = 2,
            Method = :de_rand_1_bin, TraceMode = :silent, MaxSteps = 25001)
        @test best_fitness(res) < 0.001
    end

    @testset "Fixed-dimensional problem takes precedence over search range and related params" begin
        prob = BlackBoxOptim.minimization_problem((x) -> sum(x), "no name", (10.0, 20.0), 3)
        res = bboptimize(prob; SearchRange = (0.0, 2.0), NumDimensions = 2, TraceMode = :silent)
        xbest = best_candidate(res)
        @test length(xbest) == 3
        @test xbest[1] >= 10.0
        @test xbest[2] >= 10.0
        @test xbest[3] >= 10.0
    end

    @testset "Fault Handling" begin
        #@testset "anydimensional problem" begin
        #    @test_throws ArgumentError bboptimize(BlackBoxOptim.anydim_problem("dummy", (x) -> sum(x), 0.0:1.0))
        #end

        @testset "func & SearchRange given but no NumDimensions specified" begin
            @test_throws KeyError BlackBoxOptim.setup_problem(sum, ParamsDict(:SearchRange => (0.0, 1.0)))
        end
    end

#  @testset "restarting an optimizer again and again should gradually improve" begin
#    optimizer, problem, params = BlackBoxOptim.setup_bboptimize(rosenbrock2d,
#      {:SearchRange => (-5.0, 5.0), :NumDimensions => 2,
#        :MaxSteps => 10, :TraceMode => :silent})
#    best10, fitness10, termination_reason10, elapsed_time10, params, num_evals10 = #BlackBoxOptim.run_optimizer(optimizer, problem, params);
#    best20, fitness20, termination_reason20, elapsed_time20, params, num_evals20 = #BlackBoxOptim.run_optimizer(optimizer, problem, params);
#    params[:MaxSteps] = 980
#    best1000, fitness1000, termination_reason1000, elapsed_time1000, params, num_evals1000 = #BlackBoxOptim.run_optimizer(optimizer, problem, params);
#    params[:MaxSteps] = 1000
#    fitness10000 = best10000 = elapsed_time1000b = 1 # Just so saved outside of loop body...
#    for i in 1:9
#      best10000, fitness10000, termination_reason10000, elapsed_time1000b, params, num_evals10000 = #BlackBoxOptim.run_optimizer(optimizer, problem, params);
#    end
#
#    # Fitness is not worse in sub-sequent runs
#    @test fitness10 >= fitness20 >= fitness1000 >= fitness10000
#
#    # and it should (almost) always be better after 10000 steps than after 10:
#    @test fitness10 > fitness10000
#  end
end
