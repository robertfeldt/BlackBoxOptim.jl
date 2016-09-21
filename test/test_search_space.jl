@testset "Search space" begin
    @testset "in()" begin
        for i in 1:NumTestRepetitions
            reps = rand(1:10)
            ss1 = symmetric_search_space(reps, (0.0, 1.0))
            ind = rand_individual(ss1)
            for j in 1:reps
                @test (mins(ss1)[j] <= ind[j] <= maxs(ss1)[j]) 
            end
        end
    end

    @testset "Symmetric search space with default range" begin
        ss1 = symmetric_search_space(1)
        @test numdims(ss1) == 1
        @test ranges(ss1) == [(0.0, 1.0)]
        @test range_for_dim(ss1,1) == (0.0, 1.0)

        for i in 1:NumTestRepetitions
            ind = rand_individual(ss1)
            @test size(ind) == (1,)
            @test in(ind, ss1)
        end

        ss3 = symmetric_search_space(3)
        @test numdims(ss3) == 3
        @test ranges(ss3) == [(0.0, 1.0), (0.0, 1.0), (0.0, 1.0)]
        @test range_for_dim(ss3,1) == (0.0, 1.0)
        @test range_for_dim(ss3,2) == (0.0, 1.0)
        @test range_for_dim(ss3,3) == (0.0, 1.0)

        for i in 1:NumTestRepetitions
            ind = rand_individual(ss3)
            @test size(ind) == (3,)
            @test in(ind, ss3)
        end
    end

    @testset "SymmetricSearchSpace with given range" begin
        ss1 = symmetric_search_space(1, (-1.0, 1.0))
        @test numdims(ss1) == 1
        @test ranges(ss1) == [(-1.0, 1.0)]
        @test range_for_dim(ss1,1) == (-1.0, 1.0)

        for i in 1:NumTestRepetitions
            reps = rand(1:100)
            range = (rand(), rand())
            ss = symmetric_search_space(reps, range)
            @test numdims(ss) == reps
            @test all([(dr == range) for dr in ranges(ss)])
        end
    end

    @testset "rand_individual is within the search space" begin
        for i in 1:NumTestRepetitions
            reps = rand(1:100)
            mm = sort(rand(2,1), 1)
            range = (mm[1], mm[2])
            ss = symmetric_search_space(reps, range)
            ind = rand_individual(ss)
            @test length(ind) == numdims(ss)
            @test in(ind, ss)
        end
    end

    @testset "rand_individuals creates many individuals and all are within the search space" begin
        for i in 1:NumTestRepetitions
            reps = rand(1:10)
            mm = sort(rand(2,1), 1)
            range = (mm[1], mm[2])
            ss = symmetric_search_space(reps, range)
            numinds = rand(1:10)
            inds = rand_individuals(ss, numinds)
            @test size(inds,1) == numdims(ss)
            @test size(inds,2) == numinds
            for j in 1:numinds
                @test in(inds[:,j], ss)
            end
        end
    end

    @testset "rand_individuals correctly handles column-wise generation in assymetric search spaces" begin
        for i in 1:NumTestRepetitions÷10
            numdimensions = rand(1:13)
            minbounds = rand(numdimensions)
            ds = rand(1:10, numdimensions) .* rand(numdimensions)
            maxbounds = minbounds .+ ds
            parambounds = collect(zip(minbounds, maxbounds))
            ss = RangePerDimSearchSpace(parambounds)
            @test mins(ss) == minbounds
            @test maxs(ss) == maxbounds
            @test round(deltas(ss), 6) == round(ds, 6)

            # Now generate 100 individuals and make sure they are all within bounds
            inds = rand_individuals(ss, 100)
            @test size(inds, 2) == 100
            for i in 1:size(inds, 2)
                for d in 1:numdimensions
                    @test (minbounds[d] <= inds[d,i] <= maxbounds[d])
                end
            end
        end
    end

    @testset "RangePerDimSearchSpace" begin
        ss = RangePerDimSearchSpace([(0.0, 1.0)])
        @test mins(ss) == [0.0]
        @test maxs(ss) == [1.0]
        @test deltas(ss) == [1.0]

        ss = RangePerDimSearchSpace([(0.0, 1.0), (0.5, 10.0)])
        @test mins(ss) == [0.0, 0.5]
        @test maxs(ss) == [1.0, 10.0]
        @test deltas(ss) == [1.0, 9.5]
    end

    @testset "rand_individuals_lhs samples in LHS intervals" begin
        ss = RangePerDimSearchSpace([(0.0, 1.0), (2.0, 3.0), (4.0, 5.0)])

        inds = rand_individuals_lhs(ss, 2)
        @test size(inds, 1) == 3
        @test size(inds, 2) == 2

        sorted = sort(inds, 2) # Sort per row --> in their ordered intervals
        @test (0.0 <= sorted[1,1] <= 0.5)
        @test (0.5 <= sorted[1,2] <= 1.0)

        @test (2.0 <= sorted[2,1] <= 2.5)
        @test (2.5 <= sorted[2,2] <= 3.0)

        @test (4.0 <= sorted[3,1] <= 4.5)
        @test (4.5 <= sorted[3,2] <= 5.0)
    end

    @testset "feasible finds feasible points in the search space" begin
        ss = RangePerDimSearchSpace([(0.0, 1.0), (2.0, 3.0), (4.0, 5.0)])

        # We use the double transpose below to ensure the actual and expected
        # values have the same type (matrices, not vectors).
        @test BlackBoxOptim.feasible([1.1, 2.0, 4.0], ss) == [1.0, 2.0, 4.0]
        @test BlackBoxOptim.feasible([1.1, 3.0, 4.0], ss) == [1.0, 3.0, 4.0]
        @test BlackBoxOptim.feasible([1.1, 2.0, 5.0], ss) == [1.0, 2.0, 5.0]
        @test BlackBoxOptim.feasible([1.1, 3.0, 5.0], ss) == [1.0, 3.0, 5.0]

        @test BlackBoxOptim.feasible([-0.1, 2.0, 4.0], ss) == [0.0, 2.0, 4.0]
        @test BlackBoxOptim.feasible([-0.1, 3.0, 4.0], ss) == [0.0, 3.0, 4.0]
        @test BlackBoxOptim.feasible([-0.1, 2.0, 5.0], ss) == [0.0, 2.0, 5.0]
        @test BlackBoxOptim.feasible([-0.1, 3.0, 5.0], ss) == [0.0, 3.0, 5.0]

        @test BlackBoxOptim.feasible([0.0, 1.9, 4.0], ss) == [0.0, 2.0, 4.0]
        @test BlackBoxOptim.feasible([0.0, 1.9, 4.0], ss) == [0.0, 2.0, 4.0]
        @test BlackBoxOptim.feasible([1.0, 1.9, 5.0], ss) == [1.0, 2.0, 5.0]
        @test BlackBoxOptim.feasible([1.0, 1.9, 5.0], ss) == [1.0, 2.0, 5.0]

        @test BlackBoxOptim.feasible([0.0, 3.3, 4.0], ss) == [0.0, 3.0, 4.0]
        @test BlackBoxOptim.feasible([0.0, 3.2, 4.0], ss) == [0.0, 3.0, 4.0]
        @test BlackBoxOptim.feasible([1.0, 3.1, 5.0], ss) == [1.0, 3.0, 5.0]
        @test BlackBoxOptim.feasible([1.0, 3.9, 5.0], ss) == [1.0, 3.0, 5.0]

        @test BlackBoxOptim.feasible([-0.4, 3.3, 14.5], ss) == [0.0, 3.0, 5.0]
    end

    @testset "diameters" begin
        ss = RangePerDimSearchSpace([(0.0, 1.0), (2.0, 3.0), (4.0, 5.0)])
        diams = diameters(ss)

        @test length(diams) == 3
        @test diams[:] == [1.0, 1.0, 1.0]
    end

    @testset "concat(ss1, ss2)" begin
        ss1 = RangePerDimSearchSpace([(0.0, 1.0), (2.0, 3.0), (4.0, 5.0)])
        ss2 = RangePerDimSearchSpace([(6.0, 7.0), (8.0, 9.0)])

        sscat = vcat(ss1, ss2)
        @test numdims(sscat) == 5
        @test mins(sscat) == [0.0, 2.0, 4.0, 6.0, 8.0]
        @test maxs(sscat) == [1.0, 3.0, 5.0, 7.0, 9.0]
        @test deltas(sscat) == fill(1.0, 5)
    end
end
