haltonnumber = BlackBoxOptim.Utils.haltonnumber
haltonsequence = BlackBoxOptim.Utils.haltonsequence

# The Halton sequence with base 2 according to http://en.wikipedia.org/wiki/Halton_sequence
HaltonSequence2 = [
    1/2, 1/4, 3/4, 1/8, 5/8, 3/8, 7/8, 1/16, 9/16
]

# The Halton sequence with base 2 according to http://en.wikipedia.org/wiki/Halton_sequence
HaltonSequence3 = [
    1/3, 2/3, 1/9, 4/9, 7/9, 2/9, 5/9, 8/9, 1/27
]

@testset "Halton numbers and sequence" begin

    @testset "Halton numbers" begin

        map(1:length(HaltonSequence2)) do i
            @test isapprox(haltonnumber(2, i), HaltonSequence2[i])
        end

        map(1:length(HaltonSequence3)) do i
            @test isapprox(haltonnumber(3, i), HaltonSequence3[i])
        end
    
    end

    @testset "Halton sequences" begin

        hseq2 = haltonsequence(2, length(HaltonSequence2))
        @test length(hseq2) == length(HaltonSequence2)
        map(1:length(HaltonSequence2)) do i
            @test isapprox(hseq2[i], HaltonSequence2[i])
        end

        hseq3 = haltonsequence(3, length(HaltonSequence3))
        @test length(hseq3) == length(HaltonSequence3)
        map(1:length(HaltonSequence3)) do i
            @test isapprox(hseq3[i], HaltonSequence3[i])
        end

    end

end