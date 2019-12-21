using BlackBoxOptim: general_stop_reason

@testset "general_stop_reason" begin
    problem = BlackBoxOptim.example_problems["Sphere"]

    @testset "Within fitness tolerance" begin
        # Set a low MinDeltaFitnessTolerance so it does not interfere.
        res = bboptimize(problem; NumDimensions = 1, TraceMode = :silent, 
            FitnessTolerance = 1e-1, MinDeltaFitnessTolerance = 1e-100)
        @test general_stop_reason(res) == "Within fitness tolerance of optimum"
    end

    @testset "Delta fitness below tolerance" begin
        # Set a low FitnessTolerance so it does not interfere.
        res = bboptimize(problem; NumDimensions = 1, TraceMode = :silent, 
            MinDeltaFitnessTolerance = 1e-1, FitnessTolerance = 1e-100)
        @test general_stop_reason(res) == "Delta fitness below tolerance"
    end
end