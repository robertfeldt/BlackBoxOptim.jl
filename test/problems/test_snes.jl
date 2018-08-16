include("common.jl")

@testset "Optimize single objective problems in 5, 10, 30 and 100 dimensions with sNES" begin

@testset "$problem" for problem in ["Sphere", "Schwefel2.22", "Schwefel2.22"]
    p = BlackBoxOptim.Problems.examples[problem]

    @test fitness_for_opt(p, 5, 20,  1e2, separable_nes) < 0.01
    @test fitness_for_opt(p, 5, 20,  5e2, xnes) < 0.01

    @test fitness_for_opt(p, 10, 20, 3e2, separable_nes) < 0.01
    @test fitness_for_opt(p, 10, 20, 1e3, xnes) < 0.01

    @test fitness_for_opt(p, 30, 25, 1e3, separable_nes) < 0.01
    @test fitness_for_opt(p, 30, 25, 4e3, xnes) < 0.01

    @test fitness_for_opt(p, 100, 25, 5e3, separable_nes) < 0.01
    # Cannot run xnes in 100 dimensions since it scales badly.
end

@testset "Schwefel1.2" begin
    p = BlackBoxOptim.Problems.examples["Schwefel1.2"]

    @test fitness_for_opt(p, 30, 50, 4e3, separable_nes) < 10.0
    @test fitness_for_opt(p, 30, 50, 8e3, xnes) < 10.0
end

@testset "Rosenbrock" begin
    p = BlackBoxOptim.Problems.examples["Rosenbrock"]

    @test fitness_for_opt(p, 30, 40, 2e4, separable_nes) < 100.0
    @test fitness_for_opt(p, 30, 40, 5e4, xnes) < 100.0

    @test fitness_for_opt(p, 50, 40, 3e4, separable_nes) < 100.0
end

end
