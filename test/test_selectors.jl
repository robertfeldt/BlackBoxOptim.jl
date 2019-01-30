@testset "Selection operators" begin

ss = RectSearchSpace(1, (0.0, 10.0))
fake_problem = FunctionBasedProblem(x -> 0.0, "test_problem", MinimizingFitnessScheme, ss)
fake_pop = reshape(collect(1.0:10.0), (1, 10))

@testset "SimpleSelector" begin
    @test popsize(fake_pop) == 10
    sel = BlackBoxOptim.SimpleSelector()
    for i in 1:NumTestRepetitions
        numSamples = rand(1:8)
        sampled = BlackBoxOptim.select(sel, fake_pop, numSamples)

        @test length(sampled) == numSamples

        # All sampled indices are indices into the population
        @test all(index -> in(index, 1:popsize(fake_pop)), sampled)
    end
end

@testset "RadiusLimitedSelector" begin
    sel = BlackBoxOptim.RadiusLimitedSelector(20)

    for i in 1:NumTestRepetitions
        numSamples = rand(1:sel.radius)
        sampled = BlackBoxOptim.select(sel, fake_pop, numSamples)

        @test length(sampled) == numSamples

        # All sampled indices are indices into the population
        @test all(index -> in(index, 1:popsize(fake_pop)), sampled)

        mini, maxi = minimum(sampled), maximum(sampled)
        if (maxi - mini) > max(numSamples+2, sel.radius)
            @test mini + popsize(fake_pop) <= maxi + sel.radius
        end
    end
end

@testset "Tournament" begin
    sel = BlackBoxOptim.TournamentSelector(MinimizingFitnessScheme, 3)
    fake_pop = FitPopulation(reshape(collect(1.0:10.0), (1, 10)),
                             NaN, collect(1.0:10.0))

    @testset "tournament()" begin
        for i in 1:NumTestRepetitions
            tourSize = rand(0:popsize(fake_pop))
            cand_ixs = findall(bitrand(popsize(fake_pop)))

            winner_ix = BlackBoxOptim.tournament(sel, fake_pop, cand_ixs)

            if isempty(cand_ixs)
                @test winner_ix == 0
            else
                @test in(winner_ix, 1:popsize(fake_pop))
                @test fitness(fake_pop, winner_ix) == minimum(fake_pop.fitness[cand_ixs])
            end
        end
    end

    @testset "select()" begin
        for i in 1:NumTestRepetitions
            numSamples = rand(1:3)
            sampled = BlackBoxOptim.select(sel, fake_pop, numSamples)

            @test length(sampled) == numSamples

            # All sampled indices are indices into the population
            @test all(index -> in(index, 1:popsize(fake_pop)), sampled)
        end
    end
end

end
