include("helper.jl")

const ExamplesDir = joinpath(dirname(@__FILE__()), "..", "examples")
examples = sort(filter(fn -> isfile(joinpath(ExamplesDir, fn)) && ismatch(r"\.jl$", fn), 
    readdir(ExamplesDir)))

function run_examples_as_test(examplefile::String, dirpath::String)
    testcode = quote
        print_with_color(:blue, "Running example: ", $examplefile, "\n")
        @testset $examplefile begin
            err = nothing
            try
                include(joinpath($dirpath, $examplefile))
            catch _err
                err = _err
            end
            @test err == nothing
        end
    end
    eval(testcode)
end

@testset "Examples as tests" begin
    map(efn -> run_examples_as_test(efn, ExamplesDir), examples)
end