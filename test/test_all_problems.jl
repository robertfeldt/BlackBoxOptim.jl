@testset "Problems" begin
    @testset "Can create fixed dimensional problem with search space as array" begin
        p = BlackBoxOptim.fixeddim_problem((x)->x; search_space = [(0.0, 1.0)])
        @test typeof(p) === BlackBoxOptim.FixedDimProblem
    end

    @testset "Exception if no valid search space supplied when creating fixed dimensional problem" begin
        @test_throws BlackBoxOptim.fixeddim_problem((x)->x; search_space = (0.0, 1.0))
    end

end
