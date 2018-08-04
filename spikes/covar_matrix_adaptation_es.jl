function ellipsoid(x)
  res = 0.0
  cumsum = 0.0
  for xx in x
    cumsum += xx
    res += cumsum^2
  end
  res
end

function sort_population(population, mu)
 for i=1:length(pop); fitnesses(i) = pop{i}.F; end;
 [sorted_fitnesses, index] = sort(fitnesses);
 for i=1:mu; sorted_pop{i} = pop{index(i)}; end
end

function cma_es(p)
  lambda = 12
  mu = 3


end

function cma_es(p;
  trace = true,

  # Stopping criteria related
  max_seconds = 4*numdims(p), max_evals_per_dim = 1e7,
  ftol = 1e-7, xtol = 1e-10, stol = 1e-10,
  max_rounds_without_improvement = 500,

  # Starting points, will be random unless specified
  xmean = false,

  # Algorithm specific params:
  covarMatrixSampler = CholeskyCovarSampler,
  utilitiesFunc = log_utilities,

  lambda = 12,
  mu = 3,
  decompose_covar_prob = 0.4,
  xmean = false,
  sigma = 0.05*rand(1:8)*minimum(diameters(search_space(p))),
  sigmaMin = 1e-10)

  )

  N = numdims(p)
  ss = search_space(p)
  max_evals = max_evals_per_dim * N

  tau = sqrt(N)
  tau_c = N^2
  tau_sigma = sqrt(N)

  C = covarMatrixSampler(N)
  utilities = utilitiesFunc(mu, lambda)

  s = zeros(N,1)
  s_sigma = zeros(N,1)

  xbest = xmean = rand_individual(ss)   # Current best mean value.
  fbest = eval1(xbest, p)
  num_fevals = 1
  archive = BlackBoxOptim.TopListArchive(N, 10)
  add_candidate!(archive, fbest, xbest[:], num_fevals)

  fevals_last_best = num_fevals
  next_print_covar = 100
  termination_reason = "?"

  start_time = time()

  # Now lets optimize! Ensure we run at least one iteration.
  while(true)

    if (time() - start_time) > max_seconds
      termination_reason = "Exceeded time budget"
      break
    end

    if num_fevals > max_evals
      termination_reason = "Exceeded function eval budget"
      break
    end

    # Decompose only with a given probability => saves time
    if rand() <= decompose_covar_prob
      decompose!(C)
    end

    # Generate new population
    ss = multivariate_normal_sample(C, N, lambda) # std
    z = sigma * ss
    xs = repmat(xmean, 1, lambda) + z               # n*lambda

    #
    # Evaluate fitness
    fitnesses = eval_fitnesses(p, xs, lambda)
    num_fevals += lambda

    # Check if best new fitness is best ever and print some info if tracing.
    indbest = indmin(fitnesses)
    fbest_new = fitnesses[indbest]

    if fbest_new < fbest
      xbest = xs[:, indbest]
      fbest = fbest_new
      add_candidate!(archive, fbest, xbest, num_fevals)
      fevals_last_best = num_fevals

      if fitness_is_within_ftol(p, ftol, fbest)
        termination_reason = "Within ftol"
        break
      end

      if trace
        println("$(num_fevals): Best fitness = $(fbest)")
        if num_fevals > next_print_covar
          next_print_covar = num_fevals + 100
          println("covar summary: ", sumstats(C.C, (x) -> @sprintf("%.2e", x)))
          println("sigma: ", sigma)
        end
      end
    else
      if (num_fevals - fevals_last_best) > max_rounds_without_improvement * lambda
        termination_reason = "Max rounds without improvement reached"
        break
      end
    end

    # Assign weights to the best individuals according to the utilities vector.
    weights = assign_weights(lambda, fitnesses, utilities)

    # Recombination
    rw = z * weights
    rstd =
    xmean += rw

end

function cma_recombine(w, std)
  return sum(w)/length(w), sum(std)/length(std)
end

##############################################################################
# sorts population w.r.t. the individuals fitness in ascending order
function sorted_pop = SortPop(pop, mu);
 for i=1:length(pop); fitnesses(i) = pop{i}.F; end;
 [sorted_fitnesses, index] = sort(fitnesses);
 for i=1:mu; sorted_pop{i} = pop{index(i)}; end
end

############################################################################
# initialization:
############################################################################
#Cov = Matrix{Float64}(I, n, n);           # initial covariance matrix
# initializing individual population:
#Individual.y = yInit;
#Individual.w = 0;
#Individual.std = 0;
#Individual.F = fitness(Individual.y);
#for i=1:mu; ParentPop{i} = Individual; end;
#yParent = yInit;        # initial centroid parent

###############################################################################
# evolution loop of the CMA-ES
###############################################################################
while(1)
 SqrtCov = chol(Cov)';                    # "square root" of covariance matrix
 for l = 1:lambda;                        # generate lambda offspring
  OffspringIndividual.std = randn(n,1);   # line (L1a)
  OffspringIndividual.w = sigma*(SqrtCov*OffspringIndividual.std); # line (L1a)
  OffspringIndividual.y = yParent + OffspringIndividual.w;         # line (L1b)
  OffspringIndividual.F = fitness(OffspringIndividual.y);  # determine fitness (L1c)
  OffspringPop{l} = OffspringIndividual;                   # offspring complete
 end;
 ParentPop = SortPop(OffspringPop, mu);   # sort population and take mu best
 disp(ParentPop{1}.F);                    # display best fitness in population
 Recombinant = CMArecomb(ParentPop);      # (L2) perform recombination
 yParent = yParent + Recombinant.w;       # (L2) calculate new centroid parent
 s = (1-1/tau)*s + sqrt(mu/tau*(2-1/tau))*Recombinant.w/sigma;   # line (L3)
 Cov = (1-1/tau_c)*Cov + (s/tau_c)*s';                           # line (L4)
 Cov = (Cov + Cov')/2;                    # enforce symmetry of cov matrix
 s_sigma = (1-1/tau_sigma)*s_sigma + sqrt(mu/tau_sigma*(2-1/tau_sigma))*Recombinant.std; # line (L5)
 sigma = sigma*exp((s_sigma'*s_sigma - n)/(2*n*sqrt(n)));        # line (L6)
 if (sigma < sigmaMin ) break; end;       # termination condition
end
###############################################################################
