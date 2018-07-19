if !isdefined(:TimeTestExecution)
    const TimeTestExecution = false
end

module BlackBoxOptimTests

using CSV

startclocktime = time()
include("helper.jl")

import Compat.String

my_tests = [

    "utilities/test_latin_hypercube_sampling.jl",
    "utilities/test_halton_sequence.jl",
    "utilities/test_assign_ranks.jl",

    "test_parameters.jl",
    "test_fitness.jl",
    "test_evaluator.jl",
    "test_population.jl",
    "test_bimodal_cauchy_distribution.jl",
    "test_search_space.jl",
    "test_mutation_operators.jl",
    "test_crossover_operators.jl",
    "test_selectors.jl",
    "test_embedders.jl",
    "test_frequency_adaptation.jl",
    "test_archive.jl",
    "test_epsbox_archive.jl",
    "test_optimizationresult.jl",

    "test_random_search.jl",
    "test_differential_evolution.jl",
    "test_adaptive_differential_evolution.jl",
    "test_natural_evolution_strategies.jl",

    "test_borg_moea.jl",

    "test_tracing.jl",
    "test_toplevel_bboptimize.jl",
    "test_smoketest_bboptimize.jl",

    "problems/test_single_objective.jl",

    "test_generating_set_search.jl",
    "test_direct_search_with_probabilistic_descent.jl",
]

if Main.TimeTestExecution

function get_git_remote_and_branch()
    lines = split(read(`git remote -v show`, String), "\n")
    remote = match(r"[a-z0-9]+\s+([^\s]+)", lines[1]).captures[1]
    branch = strip(read(`git rev-parse --abbrev-ref HEAD`, String))
    commit = strip(read(`git rev-parse HEAD`, String))
    return remote, branch, commit
end

gitremote, gitbranch, gitcommit = get_git_remote_and_branch()
gitstr = gitremote * "/" * gitbranch * "/" * gitcommit[1:6]
versionstr = string(VERSION)

using DataFrames

TestTimingFileName = "test/timing_testing.csv"

if isfile("test/timing_testing.csv")
    timing_data = CSV.read("test/timing_testing.csv")
else
    timing_data = DataFrame(TimeStamp = AbstractString[],
                            Julia = AbstractString[],
                            Git = AbstractString[],
                            TestFile = AbstractString[],
                            Elapsed = Float64[])
end

end

using CPUTime

numtestfiles = 0
starttime = CPUtime_us()
@testset "BlackBoxOptim test suite" begin

for t in my_tests
    if Main.TimeTestExecution
        CPUtic()

        # Including the test file runs the tests in there...
        include(t)

        elapsed = CPUtoq()
        datestr = Libc.strftime("%Y%m%d %H:%M.%S", time())
        push!(timing_data, [datestr, versionstr, gitstr, t, elapsed])
    else
        include(t)
    end
    numtestfiles += 1
    print("."); flush(STDOUT);
end
println() # So Base.Test summary is correctly aligned...
end
elapsed = float(CPUtime_us() - starttime)/1e6

if Main.TimeTestExecution
    datestr = Libc.strftime("%Y%m%d %H:%M.%S", time())
    using SHA
    hash = bytes2hex(sha512(join(map(fn -> read(open(joinpath("test", fn)), String), my_tests))))[1:16]
    push!(timing_data, [datestr, versionstr, gitstr, "TOTAL TIME for $(length(my_tests)) test files, $(hash)", elapsed])
    CSV.write(TestTimingFileName, timing_data)
    println("Wrote $(nrow(timing_data)) rows to file $TestTimingFileName")
end

elapsedclock = time() - startclocktime
println("Tested $(numtestfiles) files in $(round(elapsedclock, digits=1)) seconds.")

end # module BlackBoxOptimTests
