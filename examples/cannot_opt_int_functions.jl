using BlackBoxOptim

# The functions being optimized need to return a float so trying
# to optimize this function throws an error:
function int_returning_opt_func(x::Vector{Float64})
    round(Int, sumabs(x))
end

# Let's try it
try
    bboptimize(int_returning_opt_func;
        SearchRange = (-5.0, 5.0), NumDimensions = 3, TraceMode = :silent)
catch err
    println("Error raised: ", err)
end

# Solution is to just return a Float64 instead:
fixed_func(x::Vector{Float64}) = float(int_returning_opt_func(x))

res = bboptimize(fixed_func;
        SearchRange = (-5.0, 5.0), NumDimensions = 3)
println("Best fitness: ", best_fitness(res))