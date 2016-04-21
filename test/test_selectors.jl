facts("Selection operators") do

ss = symmetric_search_space(1, (0.0, 10.0))
fake_problem = FunctionBasedProblem(x -> 0.0, "test_problem", MinimizingFitnessScheme, ss)
fake_pop = collect(1.0:10.0)'

context("SimpleSelector") do
  @fact popsize(fake_pop) --> 10
  sel = BlackBoxOptim.SimpleSelector()
  for i in 1:NumTestRepetitions
    numSamples = rand(1:8)
    sampled = BlackBoxOptim.select(sel, fake_pop, numSamples)

    @fact length(sampled) --> numSamples

    # All sampled indices are indices into the population
    @fact all(index -> in(index, 1:popsize(fake_pop)), sampled) --> true
  end
end

context("RadiusLimitedSelector") do
  sel = BlackBoxOptim.RadiusLimitedSelector(20)

  for i in 1:NumTestRepetitions
    numSamples = rand(1:sel.radius)
    sampled = BlackBoxOptim.select(sel, fake_pop, numSamples)

    @fact length(sampled) --> numSamples

    # All sampled indices are indices into the population
    @fact all(index -> in(index, 1:popsize(fake_pop)), sampled) --> true

    mini, maxi = minimum(sampled), maximum(sampled)
    if (maxi - mini) > max(numSamples+2, sel.radius)
      @fact mini + popsize(fake_pop) --> less_than_or_equal(maxi + sel.radius)
    end
  end
end

context("Tournament") do
    sel = BlackBoxOptim.TournamentSelector(MinimizingFitnessScheme, 3)
    fake_pop = FitPopulation(collect(1.0:10.0)', NaN, collect(1.0:10.0))

    context("tournament()") do
        for i in 1:NumTestRepetitions
            tourSize = rand(0:popsize(fake_pop))
            cand_ixs = find(bitrand(popsize(fake_pop)))

            winner_ix = BlackBoxOptim.tournament(sel, fake_pop, cand_ixs)

            if isempty(cand_ixs)
                @fact winner_ix --> 0
            else
                @fact in(winner_ix, 1:popsize(fake_pop)) --> true
                @fact fitness(fake_pop, winner_ix) --> minimum(fake_pop.fitness[cand_ixs])
            end
        end
    end

    context("select()") do
        for i in 1:NumTestRepetitions
            numSamples = rand(1:3)
            sampled = BlackBoxOptim.select(sel, fake_pop, numSamples)

            @fact length(sampled) --> numSamples

            # All sampled indices are indices into the population
            @fact all(index -> in(index, 1:popsize(fake_pop)), sampled) --> true
        end
    end
end

end
