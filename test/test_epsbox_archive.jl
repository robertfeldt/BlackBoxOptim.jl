@testset "EpsBoxArchive" begin
    # check that frontier elements are mutually nondominated and
    # ther are no epsilon-index duplicates
    function check_frontier(a::EpsBoxArchive)
        for (i, el_i) in enumerate(pareto_frontier(a))
            for (j, el_j) in enumerate(pareto_frontier(a))
                if j > i
                    comp, index_match = hat_compare(fitness(el_i), fitness(el_j), fitness_scheme(a))
                    if comp != 0 || index_match
                        @warn "Dominated pair: fitness($i)=$(fitness(el_i)) fitness($j)=$(fitness(el_j))"
                        return comp, index_match # return if violating pair is found
                    end
                end
            end
        end
        return (0, false) # all non-dominated, no index matches
    end

    scheme = EpsBoxDominanceFitnessScheme{2}(0.1, is_minimizing=true)

    @testset "EpsBoxFrontierIndividual" begin
        arch_indiv = BlackBoxOptim.EpsBoxFrontierIndividual(
                        convert(IndexedTupleFitness, (2.21, 1.12), scheme), [3.0, 2.0], 3, 5, 2)
        @test arch_indiv.fitness.orig == (2.21, 1.12)
        @test arch_indiv.fitness.index == (22, 11)
        @test arch_indiv.timestamp > 0
        @test arch_indiv.tag == 3
        @test arch_indiv.num_fevals == 5
        @test arch_indiv.n_restarts == 2
    end

    @testset "Constructing a small archive and adding to it" begin
        a = EpsBoxArchive(scheme, max_size=100)

        @test capacity(a) == 100
        @test length(a)   == 0
        @test isempty(a)
        @test best_candidate(a) === nothing
        @test isnafitness(best_fitness(a), fitness_scheme(a))
        @test BlackBoxOptim.noprogress_streak(a) == 0
        @test BlackBoxOptim.tagcounts(a) == Dict{Int,Int}()

        BlackBoxOptim.add_candidate!(a, convert(IndexedTupleFitness, (2.21, 1.12), scheme), [0.0, 5.0], 1, 8)
        @test capacity(a)         == 100
        @test length(a)           == 1
        @test !isempty(a)
        @test best_fitness(a) == fitness(first(pareto_frontier(a)))
        @test best_candidate(a) === params(first(pareto_frontier(a)))
        @test BlackBoxOptim.best_front_elem(a) == first(pareto_frontier(a))
        @test first(pareto_frontier(a)).num_fevals == 8
        @test BlackBoxOptim.noprogress_streak(a) == 0
        @test BlackBoxOptim.tagcounts(a) == Dict(1=>1)

        # dominated candidate is ignored and best candidate unchanged
        BlackBoxOptim.add_candidate!(a, (3.21, 1.12), [1.0, 5.0], 2)
        @test capacity(a)         == 100
        @test length(a)           == 1
        @test BlackBoxOptim.best_front_elem(a) == first(pareto_frontier(a))
        @test first(pareto_frontier(a)).num_fevals == 8
        @test BlackBoxOptim.noprogress_streak(a) == 1
        @test BlackBoxOptim.tagcounts(a) == Dict(1=>1)

        # non-dominated candidate is appended and best candidate updated
        BlackBoxOptim.add_candidate!(a, (1.21, 2.11), [1.2, 3.0], 3)
        @test capacity(a)         == 100
        @test length(a)           == 2
        @test findfirst(a.frontier, SI.Point((12,21))) !== nothing
        candi_leaf, candi_ix = findfirst(a.frontier, SI.Point((12,21)))
        @test candi_leaf[candi_ix].num_fevals == 3
        @test best_fitness(a).orig == (1.21, 2.11)
        @test best_candidate(a) == [1.2, 3.0]
        @test BlackBoxOptim.best_front_elem(a) == candi_leaf[candi_ix]
        @test BlackBoxOptim.noprogress_streak(a) == 0
        @test BlackBoxOptim.tagcounts(a) == Dict(1=>1, 3=>1)

        # another non-dominated candidate (no best update)
        BlackBoxOptim.add_candidate!(a, (0.52, 3.15), [1.5, 3.0], 3)
        @test capacity(a)         == 100
        @test length(a)           == 3
        @test best_fitness(a).orig == (1.21, 2.11)
        @test BlackBoxOptim.noprogress_streak(a) == 0
        @test BlackBoxOptim.tagcounts(a) == Dict(1=>1, 3=>2)

        # the candidate that dominates all current frontier
        BlackBoxOptim.add_candidate!(a, (0.22, 0.22), [1.5, 3.0], 5)
        @test capacity(a)         == 100
        @test length(a)           == 1
        @test best_fitness(a).orig == (0.22, 0.22)
        @test BlackBoxOptim.noprogress_streak(a) == 0
        @test BlackBoxOptim.tagcounts(a) == Dict(5=>1)

        BlackBoxOptim.add_candidate!(a, (0.21, 0.21), [1.5, 3.0], 6)
        @test capacity(a)         == 100
        @test length(a)           == 1
        @test best_fitness(a).orig == (0.21, 0.21)
        @test BlackBoxOptim.noprogress_streak(a) == 1
        @test BlackBoxOptim.tagcounts(a) == Dict(6=>1)

        @testset "noprogress_streak(since_restart=0) is reset after restart" begin
            @test BlackBoxOptim.noprogress_streak(a, since_restart=false) == 1
            @test BlackBoxOptim.noprogress_streak(a, since_restart=true) == 1
            BlackBoxOptim.notify!(a, :restart)
            @test BlackBoxOptim.noprogress_streak(a, since_restart=false) == 1
            @test BlackBoxOptim.noprogress_streak(a, since_restart=true) == 0
            BlackBoxOptim.add_candidate!(a, (0.2, 0.2), [1.5, 3.0], 6)
            @test BlackBoxOptim.noprogress_streak(a, since_restart=false) == 2
            @test BlackBoxOptim.noprogress_streak(a, since_restart=true) == 1
        end

        @testset "updating best candidate index if old and new non-dominated" begin
            a = EpsBoxArchive(scheme, max_size=100)

            BlackBoxOptim.add_candidate!(a, convert(IndexedTupleFitness, (1.0, 0.0), scheme), [0.0, 1.0], 1)
            @test best_fitness(a).orig == (1.0, 0.0)
            @test length(a)           == 1
            BlackBoxOptim.add_candidate!(a, convert(IndexedTupleFitness, (0.6, 0.6), scheme), [0.0, 2.0], 2)
            @test best_fitness(a).orig == (1.0, 0.0)
            @test best_candidate(a) == [0.0, 1.0]
            @test length(a)           == 2
            BlackBoxOptim.add_candidate!(a, convert(IndexedTupleFitness, (0.0, 1.1), scheme), [0.0, 3.0], 3)
            @test best_fitness(a).orig == (1.0, 0.0)
            @test length(a)           == 3
            BlackBoxOptim.add_candidate!(a, convert(IndexedTupleFitness, (0.0, 0.9), scheme), [0.0, 4.0], 4)
            @test best_fitness(a).orig == (0.0, 0.9)
            @test best_candidate(a) == [0.0, 4.0]
            @test length(a)           == 3
        end

        @testset "handling dulicate elements" begin
            a = EpsBoxArchive(scheme, max_size=100)

            BlackBoxOptim.add_candidate!(a, convert(IndexedTupleFitness, (1.25, 0.0), scheme), [0.0, 1.0], 1)
            @test BlackBoxOptim.best_front_elem(a) == first(pareto_frontier(a))
            @test length(a)         == 1
            BlackBoxOptim.add_candidate!(a, convert(IndexedTupleFitness, (0.6, 0.6), scheme), [0.0, 2.0], 2)
            @test best_fitness(a).orig == (0.6, 0.6)
            @test length(a)         == 2
            BlackBoxOptim.add_candidate!(a, convert(IndexedTupleFitness, (0.0, 1.3), scheme), [0.0, 3.0], 3)
            @test best_fitness(a).orig == (0.6, 0.6)
            @test best_candidate(a) == [0.0, 2.0] # no change
            @test length(a)         == 3
            BlackBoxOptim.add_candidate!(a, convert(IndexedTupleFitness, (0.6, 0.6), scheme), [0.0, 4.0], 4)
            @test best_fitness(a).orig == (0.6, 0.6)
            @test best_candidate(a) == [0.0, 2.0] # same fitness, so no change
            @test length(a)         == 3
        end
    end

    @testset "shaffer1() Pareto frontier" begin
        schaffer1(x) = (sum(abs2, x), sum(xx -> abs2(xx - 2.0), x))
        scheme = EpsBoxDominanceFitnessScheme{2}(0.1, is_minimizing=true)
        a = EpsBoxArchive(scheme)

        for t in range(0, stop=2, length=1000)
            params = [t, t]
            BlackBoxOptim.add_candidate!(a, schaffer1(params), params)
        end

        @test length(a) == 40
    end

    @testset "archive copies the individuals" begin
        a = EpsBoxArchive(scheme, max_size=100)
        @test length(a) == 0

        indiv = [1.5, 2.0]
        BlackBoxOptim.add_candidate!(a, (1.5, 2.0), indiv)
        # equal but not identical
        @test first(a.frontier).params == indiv
        @test first(a.frontier).params !== indiv

        # modify the vector
        indiv[1] = 5.0
        # test that archived version is still the same
        @test first(a.frontier).params == [1.5, 2.0]
    end

    @testset "simulate optimization" begin
        # generate random fitness from (0..2, 0..2, 0..2)
        # and lying inside the (1.0-r1)^1.5 + (1.0-r2)^1.5 + (1.0-r3)^1.5 <= 8.0
        function rfitness()
            while true
                r = ntuple(_ -> rand(), 3)
                if sum(x -> (1 - x)^1.5, r) <= 1.0
                    return 2.0 .* r
                end
            end
        end

        scheme3 = EpsBoxDominanceFitnessScheme{3}(0.1, is_minimizing=true)
        a = EpsBoxArchive(scheme3, max_size=1000)
        indiv = [0.0, 0.0]
        for i in 1:10000
            BlackBoxOptim.add_candidate!(a, rfitness(), indiv)
            # periodically check that pareto_frontier(a) is indeed
            # epsilon-index frontier
            if mod(i, 100) == 0
                @test length(a) > 0
                @test check_frontier(a) == (0, false)
            end
        end
    end
end
