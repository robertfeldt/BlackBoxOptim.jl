require("covar_matrix_self_adaptation_es.jl")

# Calls to ensure everything has been compiled on all processors.
@everywhere cmsa_es(2, rosenbrock; max_fevals = 2, trace = false) 
@everywhere cmsa_es(2, rosenbrock; max_fevals = 2, covarMatrixSampler = CholeskyCovarSampler, trace = false) 

# Parameter study
num_runs_per_dim = 25
dims = 2.^(1:6)
# Max evals per dim
max_evals_per_dim = [1000, 10000, 100000]
#max_evals_per_dim = [1000]
samplers = [EigenCovarSampler, CholeskyCovarSampler]
utils = [linear_utilities, log_utilities]
#utils = [linear_utilities]
#mus = [0.05, 0.10, 0.20, 0.30]
mus = [1, 0.05, 0.10, 0.20]
#mus = [1]
lambdas = [(n) -> int(4 + 3*floor(log(n))), (n) -> 4*n, (n) -> 4*n*n]
problems = [rosenbrock, sphere, cigar, cigtab]
#problems = [rosenbrock]

@everywhere machine1 = "macbook_pro_retina_cth"
@everywhere machine2 = "macbook_pro_ssd_bth"
@everywhere machine = machine1

# Now create parameter combinations
Params = Any[]
for(run in 1:num_runs_per_dim)
  for(prob in problems)
    for(d in dims)
      for(maxfevals in max_evals_per_dim)
        for(uf in utils)
          for(s in samplers)
            for(l in lambdas)

              if typeof(l) == Function
                l = l(d)
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

  ts = strftime("%y%m%d,%H%M%S", time())

  tic()
  x, f, fevals = cmsa_es(d, prob; 
    max_fevals = maxfevals, covarMatrixSampler = s,
    mu = mu, lambda = l,
    utilitiesFunc = uf, trace = false)
  exectime = toq()

  return "$(ts),$(run),$(prob),$(d),$(machine),$(s),$(uf),$(mu),$(l),$(fevals),$(exectime),$(f)"
end

@everywhere function one_run_and_save(args)
  res = one_run(args)
  println(STDERR, res)
  res
end

results = pmap(one_run_and_save, Params)

fh = open(strftime("parameter_studies/%y%m%d_%H%M%S_$(machine).csv", time()), "w")
println(fh, "Date,Time,RunID,Problem,N,Machine,Sampler,Utilities,Mu,Lambda,NumFuncEvals,ExecTime,Fitness")
for(i in 1:length(results))
  println(fh, results[i])
end
close(fh)
