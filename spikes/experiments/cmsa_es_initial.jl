using BlackBoxOptim

require("../covar_matrix_self_adaptation_es.jl")
require("../experiment_framework.jl")

#jfs = BlackBoxOptim.JadeFunctionSet
#myfs = {1 => jfs[1]} #, 2 => jfs[2], 3 => jfs[3], 4 => jfs[4],
#  5 => jfs[5], 6 => jfs[6]}
#myfs = jfs

p = as_fixed_dim_problem(BlackBoxOptim.example_problems["Rosenbrock"], 2^5)
problem = BlackBoxOptim.shifted(p)
ps = {1 => problem}

#myfs = {1 => problem}
#n = 4
#ps = BlackBoxOptim.as_fixed_dim_problem_set(myfs, n)

sf(p) = begin
  n = numdims(p)
  lambda = 4 + 3 * int(ceil(log10(n)))
  mu = 2 # int(ceil(lambda/2))
  cmsa_es(p, mu = mu, lambda = lambda)
end

@time repeated_runs(sf, ps, 8; experiment = "cmsa_es_initial")

#f = sphere
#f = deceptive_cuccu2011(15, 2)
#f = griewank
#@time ts, fbs, fes = compare_params([
#  (f, 10, 4 + 3 * int(ceil(log(10))), 2),
#  (f, 10, 4 + 3 * int(ceil(log(10))), 6),
#  (f, 10, 4 * 10, 2), 
#  (f, 10, 4 * 10, 4), 
#  (f, 10, 4 * 10, 20), 
#  (f, 10, 4 * 10 * 10, 2), 
#  (f, 10, 4 * 10 * 10, 200), 
#  ],
#  ((f, n, lambda, mu) -> 
#    cmsa_es(n, xtransform(n, f); mu = mu, lambda = lambda,
#      known_fmin = 0.0, ftol = 1e-8)),
#  10
#);
