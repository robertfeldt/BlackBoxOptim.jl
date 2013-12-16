require("covar_matrix_self_adaptation_es.jl")

# Calls to ensure everything has been compiled on all processors.
@everywhere cmsa_es(2, rosenbrock; max_fevals = 2, trace = false) 
@everywhere cmsa_es(2, rosenbrock; max_fevals = 2, covarMatrixSampler = CholeskyCovarSampler, trace = false) 

function rosenbrock(x)
  n = length(x)
  sum( 100*( x[2:n] - x[1:(n-1)].^2 ).^2 + ( x[1:(n-1)] - 1 ).^2 )
end

function sphere(x)
  sum(x .^ 2)
end

function cigar(x)
  x[1]^2 + 1e6 * sum(x[2:end].^2)
end

function cigtab(x)
  x[1]^2 + 1e8 * x[end]^2 + 1e4 * sum(x[2:(end-1)].^2)
end

# Parameter study
num_runs_per_dim = 1
dims = 2.^(1:1)
# Max evals per dim
#max_evals_per_dim = [1000, 10000, 100000]
max_evals_per_dim = [1000]
samplers = [EigenCovarSampler, CholeskyCovarSampler]
#utils = [linear_utilities, log_utilities]
utils = [linear_utilities]
#mus = [0.05, 0.10, 0.20, 0.30]
#mus = [1, 0.05, 0.10, 0.20]
mus = [1]
lambdas = [(n) -> int(4 + 3*floor(log(n))), (n) -> 4*n, (n) -> 4*n*n]
#problems = [rosenbrock, sphere, cigar, cigtab]
problems = [rosenbrock]

# Now create parameter combinations
Params = Any[]
for(run in 1:num_runs_per_dim)
  for(prob in problems)
    for(d in dims)
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
                  push!(Params, (run, prob, d, maxfevals*d, uf, s, mu, l))
                end
              end

            end
          end
        end
      end
    end
  end
end

@everywhere function one_run(args)
  run, prob, d, maxfevals, uf, s, mu, l = args
  tic()
  x, f, fevals = cmsa_es(d, prob; 
    max_fevals = maxfevals, covarMatrixSampler = s,
    mu = mu, lambda = l,
    utilitiesFunc = uf, trace = false)
  time = toq()

  return "$(run),$(prob),$(d),$(s),$(uf),$(mu),$(l),$(fevals),$(time),$(f)"
end

@everywhere function one_run_and_save(args)
  res = one_run(args)
  fh = open(strftime("%y%m%d_%H%M%S_$(args[2])_$(myid()).csv", time()), "w")
  println(fh, res)
  close(fh)
  res
end

results = pmap(one_run_and_save, Params[1:2])

println("\n\nRunID,N,Sampler,Utilities,Mu,Lambda,NumFuncEvals,Time,Fitness")
for(i in 1:length(results))
  println(results[i])
end
