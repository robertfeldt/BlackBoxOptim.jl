require("src/BlackBoxOptim.jl")

require("spikes/covar_matrix_self_adaptation_es.jl")
require("spikes/experiment_framework.jl")

jfs = BlackBoxOptim.JadeFunctionSet
#myfs = {1 => jfs[1], 2 => jfs[2], 3 => jfs[3], 4 => jfs[4],
#  5 => jfs[5], 6 => jfs[6]}
myfs = jfs

n = 10
ps = BlackBoxOptim.as_fixed_dim_problem_set(myfs, n)

sf(n, f) = begin
  lambda = 4 + 3 * int(ceil(log10(n)))
  mu = int(ceil(lambda/2))
  max_fevals = int(n*1e6)
  tf = BlackBoxOptim.xshifted(n, f)
  cmsa_es(n, tf, mu = mu, lambda = lambda, max_fevals = max_fevals)
end

@time repeated_runs(sf, ps, 5)

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
