@testset "sNES" begin

function assign_weights_wrapper(candi_ixs::Vector{Int})
    candidates = BlackBoxOptim.Candidate{Float64}[BlackBoxOptim.Candidate{Float64}([0.0], i, NaN) for i in candi_ixs]
    u = BlackBoxOptim.fitness_shaping_utilities_linear(length(candi_ixs))
    BlackBoxOptim.assign_weights!(similar(u), candidates, u)
end

@testset "assign_weights!()" begin
    @testset "when indices are already ordered" begin
        u = assign_weights_wrapper([1, 2])
        @test length(u) == 2
        @test isapprox(sum(u), 0.0; atol=1E-10)
        @test u[1] == 0.5
        @test u[2] == -0.5

        u = assign_weights_wrapper([1, 2, 3])
        @test length(u) == 3
        @test isapprox(sum(u), 0.0; atol=1E-10)
        @test u[1] == 1.0 - 1/3
        @test u[2] == -1/3
        @test u[3] == -1/3

        u = assign_weights_wrapper([1, 2, 3, 4])
        @test length(u) == 4
        @test isapprox(sum(u), 0.0; atol=1E-10)
        @test isapprox(u[1], 2/3-1/4)
        @test isapprox(u[2], 1/3-1/4)
        @test u[3] == -1/4
        @test u[4] == -1/4

        u = assign_weights_wrapper([1, 2, 3, 4, 5])
        @test length(u) == 5
        @test isapprox(sum(u), 0.0; atol=1E-10)
        @test isapprox(u[1], 2/3-1/5)
        @test isapprox(u[2], 1/3-1/5)
        @test u[3] == -1/5
        @test u[4] == -1/5
        @test u[5] == -1/5

        u = assign_weights_wrapper([1, 2, 3, 4, 5, 6])
        @test length(u) == 6
        @test isapprox(sum(u), 0.0; atol=1E-10)
        @test isapprox(u[1], (3/3)/(3/3+2/3+1/3)-1/6)
        @test isapprox(u[2], (2/3)/(3/3+2/3+1/3)-1/6)
        @test isapprox(u[3], (1/3)/(3/3+2/3+1/3)-1/6)
        @test u[4] == -1/6
        @test u[5] == -1/6
        @test u[6] == -1/6
    end

    @testset "when indices are not ordered" begin
        u = assign_weights_wrapper([4, 1, 2, 3])
        @test isapprox(sum(u), 0.0; atol=1E-10)
        @test isapprox(u[4], 2/3-1/4)
        @test isapprox(u[1], 1/3-1/4)
        @test u[2] == -1/4
        @test u[3] == -1/4
    end
end

end
