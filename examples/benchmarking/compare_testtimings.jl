using DataFrames, CSV

Timings = CSV.read("test/timing_testing.csv")

LatestJulia = Timings[end, :Julia]

res = by(Timings, :Julia, 
    d1 -> by(d1, :Git, 
        d2 -> by(d2, :TestFile, 
            d3 -> mean(d3[:Elapsed]))))

# Slice based on testfiles and save with their largest time difference factor:
res2 = map(unique(res[:, :TestFile])) do testfile
    slicedf = res[(res[:,:TestFile] .== testfile),:]
    # The merged column will be called x1
    mintime, maxtime = extrema(convert(Array, slicedf[:, :x1]))
    #sorton = maxtime/mintime
    sorton = abs(maxtime - mintime)
    (testfile, slicedf, sorton)
end

function print_percent_diff(origvalue, currvalue; rev=false)
    if currvalue > origvalue
        s = "+" * string(round(100.0 * (currvalue/origvalue - 1.0), digits=1)) * "%"
        color = rev ? :red : :gren
    else
        s = "-" * string(round(100.0 * (1.0 - currvalue/origvalue), digits=1)) * "%"
        color = rev ? :green : :red
    end
    printstyled(s, color=color)
end

# Sort by largest time diff:
for (testfile, slicedf, diff) in sort(res2, by = t -> t[3], rev=false)
    sorteddf = sort(slicedf, [:x1])
    mintime = sorteddf[1, :x1]
    mindesc = join(convert(Array, sorteddf[1, 1:2]), ", ")
    println(testfile, " (fastest = ", round(mintime, digits=2), " secs for ", mindesc, "):")
    for i in 2:nrow(sorteddf)
        if sorteddf[i, :Julia] == LatestJulia
            print("  * ")
        else
            print("    ")
        end
        desc = join(convert(Array, sorteddf[i, 1:(ncol(sorteddf)-1)]), ", ")
        print(desc, ": ", round(sorteddf[i, :x1], digits=2), " (")
        print_percent_diff(mintime, sorteddf[i, :x1]; rev=true)
        println(")")
    end
    println("")
end
