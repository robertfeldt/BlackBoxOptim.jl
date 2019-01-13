@testset "Mutation operators" begin
    ss = RectSearchSpace([(-1.0, 1.0), (0.0, 100.0), (-5.0, 0.0)])

    @testset "UniformMutation" begin
        gibbs = UniformMutation(ss)

        @test search_space(gibbs) == ss
        # incorrect dimensions
        @test_throws BoundsError BlackBoxOptim.apply(gibbs, 2.0, 0, 1)
        @test_throws BoundsError BlackBoxOptim.apply(gibbs, 2.0, 4, 1)

        for i in 1:numdims(ss)
            dim_range = dimrange(ss, i)
            for _ in 1:NumTestRepetitions
                t = BlackBoxOptim.apply(gibbs, 0.0, i, 1)
                @test t >= dim_range[1]
                @test t <= dim_range[2]
            end
        end
    end

    @testset "MutationClock" begin
        mc = MutationClock(UniformMutation(ss), 0.05)

        mutations_per_param = zeros(numdims(ss))
        NumReps = 2_000
        for i in 1:NumReps
            ref_ind = rand_individual(ss)
            ind = copy(ref_ind)
            BlackBoxOptim.apply!(mc, ind, 1)
            @test in(ind, ss)
            mutations_per_param[ind .!= ref_ind] .+= 1
        end
        mut_frequencies = mutations_per_param ./ NumReps
        for p in 1:numdims(ss)
            @test mut_frequencies[p] > 0.03
            @test mut_frequencies[p] < 0.07 # roughly matches the probability
        end
    end

    @testset "FixedGeneticOperatorsMixture" begin
        mx = FixedGeneticOperatorsMixture(GeneticOperator[NoMutation(), MutationClock(UniformMutation(ss), 0.05)], [0.3, 0.7])

        ref_ind = rand_individual(ss)
        n_params_mutated = 0
        for i in 1:10000
            ind = copy(ref_ind)
            BlackBoxOptim.apply!(mx, ind, 1)
            @test in(ind, ss)
            n_params_mutated += ind != ref_ind
        end
        # the number of parameters changed should roughly match the weight of MutationClock multiplied by its mutation probability
        @test n_params_mutated/numdims(ss) > 200
        @test n_params_mutated/numdims(ss) < 500
    end

    @testset "FAGeneticOperatorsMixture" begin
        mx = FAGeneticOperatorsMixture(GeneticOperator[NoMutation(), MutationClock(UniformMutation(ss), 0.05)], pmin = 0.05, pmax = 20.0)

        ref_ind = rand_individual(ss)
        n_params_mutated = 0
        for i in 1:10000
            ind = copy(ref_ind)
            sel_op, tag = BlackBoxOptim.next(mx)
            BlackBoxOptim.apply!(sel_op, ind, 1)
            @test in(ind, ss)
            is_mutated = ind != ref_ind
            n_params_mutated += is_mutated
            is_improved = rand() < 0.25 && is_mutated
            # fake improvement in 0.25 times mutation clock is applied
            BlackBoxOptim.adjust!(mx, tag, 1, is_improved ? 1.0 : 0.0, 0.0, is_mutated)
        end
        # FA should adjust frequencies to [almost] always apply mutation clock
        # the number of parameters changed should roughly match the MutationClock mutation probability
        @test n_params_mutated/numdims(ss) > 300
        @test n_params_mutated/numdims(ss) < 700
        @test frequencies(mx)[1] < 0.1
        @test frequencies(mx)[2] > 0.9
    end

end
