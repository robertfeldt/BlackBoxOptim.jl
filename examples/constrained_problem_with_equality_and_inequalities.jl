# We optimize problem in example 3 of 
# https://towardsdatascience.com/how-to-solve-a-constraint-optimization-problem-in-r-fdf5abee197b
using BlackBoxOptim

# Original fitness function:
origfitness(x::Vector{Float64}) = x[1] * x[4] * (x[1] +x[2] + x[3]) + x[3]

Dim = 4
lowbounds  = 1.0 .* ones(Dim)
highbounds = 5.0 .* ones(Dim)
SearchRange = [(lowbounds[i], highbounds[i]) for i in 1:Dim]

# Inequality constraint:
ineq_constraint(x::Vector{Float64}) = 25 - x[1]*x[2]*x[3]*x[4] # should be <= 0.0

# Equality constraint:
eq_constraint(x::Vector{Float64}) = x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2 - 40 # should be == 0.0

# From this we can define the residual function, here using an l2 norm.
# Note that we use abs for equality and max for inequality.
subpenalties(x::Vector{Float64}) = Float64[abs(eq_constraint(x)), max(0.0, ineq_constraint(x))]
l2(vals::Vector{Float64}) = sqrt(sum(v -> v^2, vals))
residual(x::Vector{Float64}) = l2(subpenalties(x))

# To use a static penalty function we need to select a penalty constant K, which is typically large.
# However, I tried for K from 10.0 and up to 1e6 and they all give the same result... The problem
# is a very simple one though.
K = 1e3
penalized_fitness(x::Vector{Float64}) = origfitness(x) + K * residual(x)

# Let's run the optimizer. This problem seems harder to run for longer.
res = bboptimize(penalized_fitness; SearchRange, MaxTime = 20.0);

# The web page from above states that the optimum is at:
#   (1.00000000, 4.74299963, 3.82114998, 1.37940829)
# so let's check what we get with BBO:
best = best_candidate(res)
@assert isapprox(best[1], 1.0, atol = 1e-4)
@assert isapprox(best[2], 4.74299963, atol = 1e-4)
@assert isapprox(best[3], 3.82114998, atol = 1e-4)
@assert isapprox(best[4], 1.37940829, atol = 1e-4)

# This seems better than what is achieved on the web page above but we also
# run many more function evaluations...