require("covar_matrix_self_adaptation_es.jl")

# Calls to ensure everything has been compiled on all processors.
@everywhere cmsa_es(2, rosenbrock; max_fevals = 2, trace = false) 
@everywhere cmsa_es(2, rosenbrock; max_fevals = 2, covarMatrixSampler = CholeskyCovarSampler, trace = false) 

# Parameter study
num_runs = 10
ns = [2, 4, 8, 16, 32, 64]
max_evals = [1000, 10000, 100000, 1000000]
samplers = [EigenCovarSampler, CholeskyCovarSampler]
utils = [linear_utilities, log_utilities]
mus = [0.05, 0.10, 0.20, 0.30]
lambdas = [8, (n) -> int(4 + 3*floor(log(n))), (n) -> 4*n, (n) -> 4*n*n]

# Now create parameter combinations
Params = Any[]
for(run in 1:num_runs)
  for(n in ns)
    for(maxfevals in max_evals)
      for(uf in utils)
        for(s in samplers)
          for(l in lambdas)
            if typeof(l) == Function
              l = l(n)
            end

            for(mu in mus)
              if typeof(mu) != Int64
                mu = int(mu * l)
              end
              if mu < l
                if mu < 1
                  mu = 1
                end
                push!(Params, (run, n, maxfevals, uf, s, mu, l))
              end
            end

          end
        end
      end
    end
  end
end

@everywhere function one_run(args)
  run, n, maxfevals, uf, s, mu, l = args
  tic()
  x, f, fevals = cmsa_es(n, rosenbrock; 
    max_fevals = maxfevals, covarMatrixSampler = s,
    mu = mu, lambda = l,
    utilitiesFunc = uf, trace = false)
  time = toq()

  return "$(run),$(n),$(s),$(uf),$(mu),$(l),$(fevals),$(time),$(f)"
end

results = pmap(one_run, Params)

println("\n\nRunID,N,Sampler,Utilities,Mu,Lambda,NumFuncEvals,Time,Fitness")
for(i in 1:length(results))
  println(results[i])
end