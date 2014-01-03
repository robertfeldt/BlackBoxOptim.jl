using BlackBoxOptim
using JSON

require("../covar_matrix_self_adaptation_es.jl")
require("../experiment_framework.jl")

n = 8
p = as_fixed_dim_problem(BlackBoxOptim.example_problems["Rosenbrock"], n)
p = BlackBoxOptim.shifted(p)
diameter = minimum(diameters(search_space(p)))

function read_matrix_from_file(filename, nrows, ncols)
  fh = open(filename, "r")
  fromR = JSON.parse(readall(fh))
  close(fh)
  reshape(fromR, nrows, ncols)
end

function write_array_to_R(ys)
  join(["c(", join(ys, ","), ")"])
end

params = read_matrix_from_file("params.txt", 7, 6)

function run_cmsa_based_on_params_in_file(filename, numrows)

  params = read_matrix_from_file("params2.txt", numrows, 6)

  ys = zeros(size(params,1), 3)

  for(i in 1:size(params,1))
    ps = params[i,:]
  
    lambda = 4 * [2, n, n*n][int(ceil(ps[1]))]
    mu = int(lambda / [2, 4, 8][int(ceil(ps[2]))])
    sampler = [EigenCovarSampler, CholeskyCovarSampler][int(ceil(ps[3]))]
    utilsFunc = [log_utilities, linear_utilities][int(ceil(ps[4]))]
    covar_learning_rate = ps[5]
    sigma = diameter * 10^ps[6]

    tic()
    xb, fb, nf, r, a = cmsa_es(p; max_evals_per_dim = 1e7,
      lambda = lambda,
      mu = mu,
      covarMatrixSampler = sampler,
      utilitiesFunc = utilsFunc,
      covar_learning_rate = covar_learning_rate,
      sigma = sigma)
    t = toq()
    ys[i,:] = [fb, nf, t]
  end
  write_array_to_R(ys[:,1])
end

sf(p) = begin
  p = BlackBoxOptim.shifted(p)
  n = numdims(p)
  #lambda = 4 + 3 * int(ceil(log(n)))
  #lambda = 4*n
  lambda = 4*n*n
  #lambda = rand( (4 + 3 * int(ceil(log(n)))):(4*n*n) )
  mu = int(ceil(lambda/4))
  #max_rounds_without_improvement = 1000
  xb, fb, nf, r, a = cmsa_es(p, lambda = lambda, mu = mu, utilitiesFunc = linear_utilities)
  return "$(lambda),$(mu),Linear", "pLambda,pMu,pUtilFunc", xb, fb, nf, r, a
end

ps = BlackBoxOptim.as_fixed_dim_problem_set(BlackBoxOptim.example_problems, 2.^[2,3,4])

@time repeated_runs(sf, ps, 5; experiment = "cmsa_es_std")
