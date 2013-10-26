module Problems

export OptimizationProblem

immutable OptimizationProblem
  name::ASCIIString

  # Objective functions
  funcs::Vector{Function}

  # Range per dimension
  range_per_dimension::(Float64, Float64)
end

examples = Dict{ASCIIString, OptimizationProblem}()

include("single_objective.jl")

end # module Problems