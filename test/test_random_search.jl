@testset "Random search" begin

@testset "ask()" begin
    for i in 1:NumTestRepetitions
        dims = rand(1:20)
        min = rand(1:123)
        range = (min * rand(), min + rand() * min)

        ss = RectSearchSpace(dims, range)
        opt = BlackBoxOptim.random_search(ss)

        res1 = BlackBoxOptim.ask(opt)

        @test length(res1) == 1
        candidate = res1[1]
        @test length(candidate.params) == dims
        @test in(candidate.params, ss)

        # Fake fitness
        candidate.fitness = 1.0
        better = BlackBoxOptim.tell!(opt, [candidate])

        @test better == 1
        @test opt.best == candidate.params
        @test opt.best_fitness == 1.0

        # Get one more and fake that it has better fitness.
        res2 = BlackBoxOptim.ask(opt)
        @test length(res2) == 1
        candidate2 = res2[1]
        @test length(candidate2.params) == dims
        @test in(candidate2.params, ss)
        candidate2.fitness = 0.5
        better2 = BlackBoxOptim.tell!(opt, [candidate2])
    end
end

end
