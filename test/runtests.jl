if !@isdefined(TimeTestExecution)
    const TimeTestExecution = false
end

module BlackBoxOptimTests

using LinearAlgebra, Random
using Printf: @printf, @sprintf
using SpatialIndexing

const SI = SpatialIndexing

TestDir = first(splitdir(@__FILE__()))

# If two arguments the second one if filename of a testset file
# listing the testfiles to use.
if length(ARGS) == 2 && isfile(ARGS[2])
    test_file_list = ARGS[2]
else
    test_file_list = joinpath(TestDir, "testset_normal.txt")
end

TestFiles = filter(fn -> isfile(joinpath(TestDir, fn)), 
    readlines(test_file_list))

function latest_changed_file(files, dir = "")
    files[first(sortperm(map(fn -> mtime(joinpath(dir, fn)), files), rev = true))]
end

# If a first argument is given it must be either:
#  all =>           run all test files in the testset
#  latestchanged => run only the latest changed file in the testset
if length(ARGS) > 0
    if ARGS[1] == "all"
        TestFiles = TestFiles # Change nothing so run them all
    elseif ARGS[1] == "latestchanged"
        TestFiles = AbstractString[latest_changed_file(TestFiles, TestDir)]
        println("Testing files: $(TestFiles)")
    end
end

startclocktime = time()
include("helper.jl")

import Compat.String

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

using CSV, DataFrames

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

starttime = CPUtime_us()
@testset "BlackBoxOptim test suite" begin

for t in TestFiles
    Main.TimeTestExecution && CPUtic()

    # Including the test file runs the tests in there...
    include(t)

    if Main.TimeTestExecution
        elapsed = CPUtoq()
        datestr = Libc.strftime("%Y%m%d %H:%M.%S", time())
        push!(timing_data, [datestr, versionstr, gitstr, t, elapsed])
    end
    print("."); flush(stdout);
end
println() # So Base.Test summary is correctly aligned...
end
elapsed = float(CPUtime_us() - starttime)/1e6

if Main.TimeTestExecution
    datestr = Libc.strftime("%Y%m%d %H:%M.%S", time())
    using SHA
    hash = bytes2hex(sha512(join(map(fn -> read(open(joinpath("test", fn)), String), TestFiles))))[1:16]
    push!(timing_data, [datestr, versionstr, gitstr, "TOTAL TIME for $(length(TestFiles)) test files, $(hash)", elapsed])
    CSV.write(TestTimingFileName, timing_data)
    println("Wrote $(nrow(timing_data)) rows to file $TestTimingFileName")
end

elapsedclock = time() - startclocktime
println("Tested $(length(TestFiles)) files in $(round(elapsedclock, digits=1)) seconds.")

end # module BlackBoxOptimTests
