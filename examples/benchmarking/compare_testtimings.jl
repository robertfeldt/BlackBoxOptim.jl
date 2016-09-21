using DataFrames

Timings = readtable("test/timing_testing.csv")

res = by(Timings, :Julia, d1 -> by(d1, :Git, d2 -> by(d2, :TestFile, d3 -> mean(d3[:Elapsed]))))

# Slice based on testfiles and save with their largest time difference factor:
res2 = map(unique(res[:, :TestFile])) do testfile
    slicedf = res[(res[:,:TestFile] .== testfile),:]
    # The merged column will be called x1
    mintime, maxtime = extrema(convert(Array, slicedf[:, :x1]))
    #sorton = maxtime/mintime
    sorton = abs(maxtime - mintime)
    (testfile, slicedf, sorton)
end

function percent_diff(origvalue, currvalue)
    if currvalue > origvalue
        "+" * string(round(100.0 * (currvalue/origvalue - 1.0), 1)) * "%"
    else
        "-" * string(round(100.0 * (1.0 - currvalue/origvalue), 1)) * "%"
    end
end

# Sort by largest time diff:
for (testfile, slicedf, diff) in sort(res2, by = t -> t[3], rev=true)
    sorteddf = sort(slicedf, cols=[:x1])
    mintime = sorteddf[1, :x1]
    mindesc = join(convert(Array, sorteddf[1, 1:2]), ", ")
    println(testfile, " (fastest = ", round(mintime, 2), " secs for ", mindesc, "):")
    for i in 2:nrow(sorteddf)
        desc = join(convert(Array, sorteddf[i, 1:(ncol(sorteddf)-1)]), ", ")
        println("  ", desc, ": ", round(sorteddf[i, :x1], 2), " (", percent_diff(mintime, sorteddf[i, :x1]), ")")
    end
    println("")
end
