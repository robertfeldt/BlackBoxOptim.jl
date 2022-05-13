# We optimize problem in example 2 of 
# https://towardsdatascience.com/how-to-solve-a-constraint-optimization-problem-in-r-fdf5abee197b
using BlackBoxOptim

# Original fitness function:
origfitness(x::Vector{Float64}) = sqrt(x[2])

# where x[2] >= 0.0. The other dimensions and upper bounds are unclear but let's go for 2 and 
# some arbitrary high bound.
Dim = 2
lowbounds = 0.0 .* ones(Dim)
HB = 1e6 # This should be Inf but let's be practical...
highbounds = HB .* ones(Dim)
SearchRange = [(lowbounds[i], highbounds[i]) for i in 1:Dim]

# Two inequality constraints:
a1, b1, a2, b2 = 2, 0, -1, 1
ineq_constraint1(x::Vector{Float64}) = (a1*x[1] + b1)^3 - x[2] # should be <= 0.0
ineq_constraint2(x::Vector{Float64}) = (a2*x[1] + b2)^3 - x[2] # should be <= 0.0

# From this we can define the residual function, here using an l2 norm:
const InEqFuncs = Function[ineq_constraint1, ineq_constraint2]
ineq_penalties(x::Vector{Float64}) = Float64[max(0.0, fn(x)) for fn in InEqFuncs] # Use max since we have <= 0.0 inequalities, change to min otherwise!
l2(vals::Vector{Float64}) = sqrt(sum(v -> v^2, vals))
residual(x::Vector{Float64}) = l2(ineq_penalties(x))

# To use a static penalty function we need to select a penalty constant K, which is typically large.
# However, I tried for K from 10.0 and up to 1e6 and they all give the same result... The problem
# is a very simple one though.
K = 1e3
penalized_fitness(x::Vector{Float64}) = origfitness(x) + K * residual(x)

# Let's run the optimizer. Problem is small/simple so we don't need to run for long, at least on my machine...
res = bboptimize(penalized_fitness; SearchRange, MaxTime = 1.0);

# The nloptr package in the web page from above got the result:
# Optimal value of objective function:  0.544331053951819
# Optimal value of controls: 0.3333333 0.2962963
# so let's check what we get with BBO:
@assert isapprox(best_fitness(res), 0.54433105, atol=1e-6)
best = best_candidate(res)
@assert isapprox(best[1], 0.3333333, atol = 1e-6)
@assert isapprox(best[2], 0.2962963, atol = 1e-6)