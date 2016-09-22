@testset "Bimodal Cauchy Distributions" begin

    @testset "sample bimodal cauchy with truncation above 1" begin
        bc = BlackBoxOptim.BimodalCauchy(0.65, 0.1, 1.0, 0.1, clampBelow0=false)

        n_gt0 = 0
        n_le1 = 0
        for _ in 1:10000
            v = rand(bc)
            (v > 0.0) && (n_gt0 += 1) # Very unlikely to be 0.0 so should be ok...
            (v <= 1.0) && (n_le1 += 1)
        end
        @test n_gt0 == 10000
        @test n_le1 == 10000
    end

    @testset "sample bimodal cauchy with truncation below 0 and above 1" begin
        bc = BlackBoxOptim.BimodalCauchy(0.1, 0.1, 0.95, 0.1)

        n_ge0 = 0
        n_le1 = 0
        for _ in 1:10000
            v = rand(bc)
            (v >= 0.0) && (n_ge0 += 1)
            (v <= 1.0) && (n_le1 += 1)
        end
        @test n_ge0 == 10000
        @test n_le1 == 10000
    end

end
