@testset "TopListArchive" begin
    @testset "TopListIndividual" begin
        # test equality
        @test isequal(BlackBoxOptim.TopListIndividual([3.0, 2.0], 2.0, 0), BlackBoxOptim.TopListIndividual([3.0, 2.0], 2.0, 0))
        @test (BlackBoxOptim.TopListIndividual([3.0, 2.0], 2.0, 0) == BlackBoxOptim.TopListIndividual([3.0, 2.0], 2.0, 0))
        @test (BlackBoxOptim.TopListIndividual([3.0, 2.0], 2.0, 0) != BlackBoxOptim.TopListIndividual([3.0, 2.0], 2.0, 0)) == false
        @test (BlackBoxOptim.TopListIndividual([3.0, 2.0], 1.0, 0) != BlackBoxOptim.TopListIndividual([3.0, 2.0], 2.0, 0))
        @test (BlackBoxOptim.TopListIndividual([1.0, 2.0], 2.0, 0) != BlackBoxOptim.TopListIndividual([3.0, 2.0], 2.0, 0))
    end

    @testset "Constructing a small archive and adding to it" begin

        a = TopListArchive(MinimizingFitnessScheme, 1, 3)

        @test capacity(a) == 3
        @test length(a)   == 0
        @test isnan(best_fitness(a))

        BlackBoxOptim.add_candidate!(a, 1.0, [0.0])
        @test capacity(a)         == 3
        @test length(a)           == 1
        @test best_fitness(a)     == 1.0
        @test best_candidate(a)   == [0.0]
        @test last_top_fitness(a) == 1.0
        @test delta_fitness(a)    == Inf

        BlackBoxOptim.add_candidate!(a, 2.0, [1.0])
        @test capacity(a)       == 3
        @test length(a)         == 2
        @test best_fitness(a)   == 1.0
        @test best_candidate(a) == [0.0]
        @test last_top_fitness(a) == 2.0
        @test delta_fitness(a)    == Inf

        BlackBoxOptim.add_candidate!(a, 0.5, [2.0])
        @test capacity(a)         == 3
        @test length(a)           == 3
        @test best_fitness(a)   == 0.5
        @test best_candidate(a) == [2.0]
        @test last_top_fitness(a) == 2.0
        @test delta_fitness(a)    == 0.5

        BlackBoxOptim.add_candidate!(a, 0.8, [4.0])
        @test capacity(a)         == 3
        @test length(a)           == 3
        @test best_fitness(a)   == 0.5
        @test best_candidate(a) == [2.0]
        @test last_top_fitness(a) == 1.0
        @test delta_fitness(a)    == 0.5

        expected = ((0.8 - 0.5) / ((1 - 0.05)^(-2/1) - 1))
        @test width_of_confidence_interval(a, 0.05) ≈ expected
        @test fitness_improvement_potential(a, 0.05) ≈ (expected / 0.50)

        BlackBoxOptim.add_candidate!(a, 0.4, [1.9])
        @test capacity(a)       == 3
        @test length(a)         == 3
        @test best_fitness(a)   == 0.4
        @test best_candidate(a) == [1.9]
        @test last_top_fitness(a) == 0.8
        @test isapprox(delta_fitness(a), 0.1)

        expected = ((0.5 - 0.4) / ((1 - 0.05)^(-2/1) - 1))
        @test width_of_confidence_interval(a, 0.05) ≈ expected
        @test fitness_improvement_potential(a, 0.05) ≈ (expected / 0.40)

        # identical candidate is not inserted
        BlackBoxOptim.add_candidate!(a, 0.5, [2.0])
        @test capacity(a)       == 3
        @test length(a)         == 3
        @test last_top_fitness(a)   == 0.8

        # different candidates with the same fitness are inserted
        BlackBoxOptim.add_candidate!(a, 0.5, [3.0])
        @test length(a)             == 3
        @test last_top_fitness(a)   == 0.5
    end

    @testset "magnitude_class for positive fitness values" begin

        @test BlackBoxOptim.magnitude_class(1.0)      == (1.0, 0.0)
        @test BlackBoxOptim.magnitude_class(10.0)     == (1.0, 1.0)
        @test BlackBoxOptim.magnitude_class(1e2)      == (1.0, 2.0)
        @test BlackBoxOptim.magnitude_class(1e102)    == (1.0, 102.0)
        @test BlackBoxOptim.magnitude_class(1e-1)     == (1.0, -1.0)
        @test BlackBoxOptim.magnitude_class(1e-2)     == (1.0, -2.0)
        @test BlackBoxOptim.magnitude_class(1e-9)     == (1.0, -9.0)
        @test BlackBoxOptim.magnitude_class(1e-12)    == (1.0, -12.0)

        @test BlackBoxOptim.magnitude_class(2.0)      == (1.0, 0.3)
        @test BlackBoxOptim.magnitude_class(20.0)     == (1.0, 1.3)
        @test BlackBoxOptim.magnitude_class(200.0)    == (1.0, 2.3)

        @test BlackBoxOptim.magnitude_class(5.0)      == (1.0, 0.6)
        @test BlackBoxOptim.magnitude_class(50.0)     == (1.0, 1.6)
        @test BlackBoxOptim.magnitude_class(500.0)    == (1.0, 2.6)

    end

    @testset "magnitude_class for negative fitness values" begin

        @test BlackBoxOptim.magnitude_class(-1.0)      == (-1.0, 0.0)
        @test BlackBoxOptim.magnitude_class(-10.0)     == (-1.0, 1.0)
        @test BlackBoxOptim.magnitude_class(-1e2)      == (-1.0, 2.0)
        @test BlackBoxOptim.magnitude_class(-1e102)    == (-1.0, 102.0)
        @test BlackBoxOptim.magnitude_class(-1e-1)     == (-1.0, -1.0)
        @test BlackBoxOptim.magnitude_class(-1e-2)     == (-1.0, -2.0)
        @test BlackBoxOptim.magnitude_class(-1e-9)     == (-1.0, -9.0)
        @test BlackBoxOptim.magnitude_class(-1e-12)    == (-1.0, -12.0)

        @test BlackBoxOptim.magnitude_class(-2.0)      == (-1.0, 0.3)
        @test BlackBoxOptim.magnitude_class(-20.0)     == (-1.0, 1.3)
        @test BlackBoxOptim.magnitude_class(-200.0)    == (-1.0, 2.3)

        @test BlackBoxOptim.magnitude_class(-5.0)      == (-1.0, 0.6)
        @test BlackBoxOptim.magnitude_class(-50.0)     == (-1.0, 1.6)
        @test BlackBoxOptim.magnitude_class(-500.0)    == (-1.0, 2.6)

    end

    @testset "archive copies the individuals" begin
        a = TopListArchive(MinimizingFitnessScheme, 2, 3)

        indiv = [0.0, 2.0]
        BlackBoxOptim.add_candidate!(a, 1.0, indiv)
        # equal but not identical
        @test (best_candidate(a) == indiv)
        @test (best_candidate(a) === indiv) == false

        # modify the vector
        indiv[1] = 5.0
        # test that archived version is still the same
        @test best_candidate(a) == [0.0, 2.0]
    end

    @testset "for maximizing fitness" begin
        a = TopListArchive(MaximizingFitnessScheme, 1, 3)

        BlackBoxOptim.add_candidate!(a, 1.0, [0.0])
        @test best_fitness(a)     == 1.0
        @test best_candidate(a)   == [0.0]
        @test last_top_fitness(a) == 1.0
        @test delta_fitness(a)    == Inf

        BlackBoxOptim.add_candidate!(a, 2.0, [1.0])
        @test best_fitness(a)   == 2.0
        @test best_candidate(a) == [1.0]
        @test last_top_fitness(a) == 1.0
        @test delta_fitness(a)    == 1.0

        BlackBoxOptim.add_candidate!(a, 0.5, [2.0])
        @test best_fitness(a)   == 2.0
        @test best_candidate(a) == [1.0]
        @test last_top_fitness(a) == 0.5
        @test delta_fitness(a)    == 1.0

        BlackBoxOptim.add_candidate!(a, 1.5, [4.0])
        @test best_fitness(a)   == 2.0
        @test best_candidate(a) == [1.0]
        @test last_top_fitness(a) == 1.0
        @test delta_fitness(a)    == 1.0

        BlackBoxOptim.add_candidate!(a, 2.5, [5.0])
        @test best_fitness(a)   == 2.5
        @test best_candidate(a) == [5.0]
        @test last_top_fitness(a) == 1.5
        @test delta_fitness(a)    == 0.5
    end
end
