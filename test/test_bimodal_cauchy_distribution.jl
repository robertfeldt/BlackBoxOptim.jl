@testset "Bimodal Cauchy Distributions" begin

    @testset "sample bimodal cauchy with truncation on one side" begin
        bc = BlackBoxOptim.bimodal_cauchy(0.65, 0.1, 1.0, 0.1)

        for i in 1:100
            v = BlackBoxOptim.sample_bimodal_cauchy(bc; truncateAbove1 = true, truncateBelow0 = false)

            @test v <= 1.0
            @test v > 0.0 # Very unlikely to be 0.0 so should be ok...
        end
    end

    @testset "sample bimodal cauchy with truncation on both sides" begin
        bc = BlackBoxOptim.bimodal_cauchy(0.1, 0.1, 0.95, 0.1)

        for i in 1:100
            v = BlackBoxOptim.sample_bimodal_cauchy(bc; truncateAbove1 = true, truncateBelow0 = true)

            @test v <= 1.0
            @test v >= 0.0
        end
    end

end
