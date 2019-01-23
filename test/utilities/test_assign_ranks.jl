@testset "Assign ranks within tolerance" begin

    @testset "Ranks correctly if none are within tolerance of each other" begin
        res = BlackBoxOptim.Utils.assign_ranks_within_tolerance([1.0, 2.0, 3.0])
        @test length(res) == 3
        @test res[1] == (1, 1.0, 1.0)
        @test res[2] == (2, 2.0, 2.0)
        @test res[3] == (3, 3.0, 3.0)

        res = BlackBoxOptim.Utils.assign_ranks_within_tolerance([3.0])
        @test length(res) == 1
        @test res[1] == (1, 3.0, 3.0)
    end

    @testset "Ranks correctly if some are within tolerance of each other" begin
        res = BlackBoxOptim.Utils.assign_ranks_within_tolerance([1.0, 1.0+1e-6, 200.0]; tolerance = 1e-5)
        @test length(res) == 3
        @test res[1] == (1, 1.0, 1.0)
        @test res[2] == (1, 1.0+1e-6, 1.0+1e-6)
        @test res[3] == (3, 200.0, 200.0)
    end

    @testset "Ranks correctly if all are within tolerance of each other" begin
        res = BlackBoxOptim.Utils.assign_ranks_within_tolerance([1.0, 1.0+1e-6, 1.0-1e-6]; tolerance = 1e-5)
        @test length(res) == 3
        @test res[1] == (1, 1.0-1e-6, 1.0-1e-6)
        @test res[2] == (1, 1.0, 1.0)
        @test res[3] == (1, 1.0+1e-6, 1.0+1e-6)
    end

    @testset "Ranks in reverse if none are within tolerance of each other" begin
        res = BlackBoxOptim.Utils.assign_ranks_within_tolerance([1.0, 20.0, 13.0]; rev = true)
        @test length(res) == 3
        @test res[1] == (1, 20.0, 20.0)
        @test res[2] == (2, 13.0, 13.0)
        @test res[3] == (3, 1.0, 1.0)
    end

    @testset "Ranks correctly for a complex example testing many aspects" begin
        values = shuffle([(-11.0, :a), (1.0, :b), (1.0+1e-6, :c), (1.0-1e-6, :d), 
            (2000, :e), (3000, :f)])
        res = BlackBoxOptim.Utils.assign_ranks_within_tolerance(values;
            tolerance = 1e-5, by = (p) -> p[1], rev = true)
        @test length(res) == length(values)
        @test res[1] == (1, (3000, :f), 3000)
        @test res[2] == (2, (2000, :e), 2000)
        @test res[3] == (3, (1.0+1e-6, :c), 1.0+1e-6)
        @test res[4] == (3, (1.0, :b), 1.0)
        @test res[5] == (3, (1.0-1e-6, :d), 1.0-1e-6)
        @test res[6] == (6, (-11.0, :a), -11.0)
    end

end
