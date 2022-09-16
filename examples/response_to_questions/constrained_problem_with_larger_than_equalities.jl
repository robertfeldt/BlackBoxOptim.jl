# This is in response to the question in
#  https://github.com/robertfeldt/BlackBoxOptim.jl/issues/205
using BlackBoxOptim

# Original fitness function (was unknown so let's just use a "random" one):
origfitness(x::Vector{Float64}) = 1.7*sqrt(x[1]) - 0.45*sqrt(x[2]*x[3]^2) - x[3]*x[1]

# Distance from fulfilling constraint 1: x[1] > x[2]
# Assign a min penalty for constraint violation and then larger penalty the further away from fulfillment
MinPenalty = 1.0
inequality1(x::Vector{Float64}) = (x[1] > x[2]) ? 0.0 : (MinPenalty + x[2] - x[1])

# Distance from fulfilling constraint 1: x[2] > x[3]
inequality2(x::Vector{Float64}) = (x[2] > x[3]) ? 0.0 : (MinPenalty + x[3] - x[2])

# where x[2] >= 0.0. The other dimensions and upper bounds are unclear but let's go for 2 and 
# some arbitrary high bound.
lowbounds =  [0.0, 0.0, 0.0]
highbounds = [1.0, 1.0, 1.0]
Dim = length(lowbounds)
SearchRange = [(lowbounds[i], highbounds[i]) for i in 1:Dim]

# From this we can define the residual function, here using an l2 norm:
const InEqFuncs = Function[inequality1, inequality2]
ineq_penalties(x::Vector{Float64}) = Float64[fn(x) for fn in InEqFuncs]
l2(vals::Vector{Float64}) = sqrt(sum(v -> v^2, vals))
residual(x::Vector{Float64}) = l2(ineq_penalties(x))

# To use a static penalty function we need to select a penalty constant K, which is typically large.
# However, I tried for K from 10.0 and up to 1e6 and they all give the same result... The problem
# is a very simple one though.
K = 1e3
penalized_fitness(x::Vector{Float64}) = origfitness(x) + K * residual(x)

# Let's run the optimizer. Problem is small/simple so we don't need to run for long, at least on my machine...
res = bboptimize(penalized_fitness; SearchRange, MaxTime = 1.0);

# Let's check constraints are fullfilled
bestsolution = best_candidate(res)
@assert inequality1(bestsolution) == 0.0
@assert bestsolution[1] > bestsolution[2]
@assert inequality2(bestsolution) == 0.0
@assert bestsolution[2] > bestsolution[3]
