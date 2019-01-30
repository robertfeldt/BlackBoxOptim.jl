@testset "Embedding operators" begin

@testset "RandomBound" begin
    @testset "does nothing if within bounds" begin
        @test apply!(RandomBound(RectSearchSpace(1, (0.0, 1.0))), [0.0], [0.0]) == [0.0]

        @test apply!(RandomBound(RectSearchSpace([(0.0, 1.0), (10.0, 15.0)])),
                     [0.0, 11.4], [0.1, 12.3] ) == [0.0, 11.4]
    end

    @testset "bounds if lower than min bound" begin
        @test apply!(RandomBound(RectSearchSpace([(0.0, 1.0)])),
                     [-0.1], [0.0]) == [0.0]

        res = apply!(RandomBound(RectSearchSpace([(0.0, 1.0)])),
                     [-0.1], [0.5])
        @test 0.0 <= res[1] <= 0.5

        res = apply!(RandomBound(RectSearchSpace([(-11.0, 1.0), (0.0, 1.0)])),
                     [-11.1, 0.5], [-10.8, 0.5])
        @test -10.8 >= res[1] >= -11.0
        @test res[2] == 0.5

        res = apply!(RandomBound(RectSearchSpace([(30.0, 60.0), (-102.0, -1.0)])),
                     [50.4, -103.1], [49.6, -101.4])
        @test res[1] == 50.4
        @test -101.4 >= res[2] >= -102.0
    end

    @testset "bounds if higher than max bound" begin
        @test apply!(RandomBound(RectSearchSpace([(0.0, 1.0)])),
                     [1.1], [1.0]) == [1.0]

        res = apply!(RandomBound(RectSearchSpace([(-10.0, 96.0)])), [97.0], [95.0])
        @test 95.0 <= res[1] <= 96.0
    end

    @testset "MixedPrecisionRectSearchSpace" begin
        ss = RectSearchSpace([(0.0, 1.0), (1.0, 2.0), (2.0, 3.0)], dimdigits=[-1, 0, 2])

        @test apply!(RandomBound(ss), [0.0, 1.4, 2.1], [0.0, 1.45, 2.0]) == [0.0, 1.0, 2.1]
        @test apply!(RandomBound(ss), [0.0, 5.0, 0.123], [0.0, 1.5, 2.0]) == [0.0, 2.0, 2.0]
        @test apply!(RandomBound(ss), [1.12345, 1.12345, 2.1234], [1.0, 2.0, 3.0]) == [1.0, 1.0, 2.12]
    end
end

end
