# Idea: Use LASSO for screening the active factors in an optimization
# problem.

using GLMNet

# Say that we want to optimize a function which has many factors but where
# only very few of them are active in calculating the outcome.
N = 1000
K = 3
active_factors = shuffle(collect(1:N))[1:K]
beta_intercept = 10.0
beta_slope = 1.0
betas = beta_intercept + beta_slope * randn(K, 1)
ofunc(x) = sum(x[active_factors] .* betas)

# Now lets sample and calc the function values for samples
S = int(10*K*log10(N));
xs = randn(N, S);
y = zeros(S)
[(y[i] = ofunc(xs[:,i])) for i in 1:S]

# Now do lasso regression:
#path = glmnet(xs', y; standardize = false)
# cv = glmnetcv(xs', y; standardize = false)

function active_factors_from_lasso(xs, y; alpha = 1.0)
  n = size(xs, 2)
  # Find the best lambda using cross-validation. We limit the size of the predictors.
  #cv = glmnetcv(xs, y) # ; constraints = hcat(-1000.0*ones(n), 1000.0*ones(n)))
  cv = glmnetcv(xs, y; standardize = false, alpha = alpha)
  num_steps = size(cv.path.betas, 2)
  res = Any[]
  for(i in 1:num_steps)
    pbetas = cv.path.betas[:, i]
    pbwi = [(pbetas[i], i) for i in 1:length(pbetas)]
    push!(res, filter((p) -> abs(p[1]) > 0.001, pbwi))
  end
  return res
end

# Count how many times each factor was found by the lasso.
function most_active_factors_from_lasso(xs, y; alpha = 1.0)
  lasso_factors = active_factors_from_lasso(xs, y; alpha = alpha)
  factor_counts = Dict{Int64, Int64}()
  for(i in 1:length(lasso_factors))
    factors = map((t) -> t[2], lasso_factors[i])
    for(f in factors)
      if haskey(factor_counts, f)
        factor_counts[f] += 1
      else
        factor_counts[f] = 1
      end
    end
  end
  return factor_counts, lasso_factors
end

function eval_if_lasso_finds_active_factors(N, K, S, reps = 100; 
  beta_intercept = 0.0, beta_slope = 1.0, alpha = 1.0)

  active_factors = shuffle(collect(1:N))[1:K]
  print("Active factors:"); show(active_factors)
  betas = beta_intercept + beta_slope * randn(K, 1)
  print("\nBetas:"); show(betas)

  objective_func(x) = begin
    active = x[active_factors]
    sum(active .* betas)
  end

  # We want to count how many times the predicted factors (ranked based on how 
  # often they were in the lasso factor set) are among the actually active factors.
  factor_found_counts = Dict{Int64, Int64}()

  for(rep in 1:reps)
    # Now lets sample
    xs = randn(N, S)

    # And calc the function values for them:
    y = zeros(S)
    for(i in 1:S)
      y[i] = objective_func(xs[:,i])
    end

    # Get the active factors from lasso reg
    factor_counts, factor_path = most_active_factors_from_lasso(xs', y; alpha = alpha)
    print("\nFactor counts:"); show(factor_counts)

    # Rank factors based on how often they were predicted as active
    factors_ranked = map((t) -> t[1], sort(collect(factor_counts), by = (t) -> t[2], rev = true))

    for(i in 1:length(factors_ranked))
      if in(factors_ranked[i], active_factors)
        if haskey(factor_found_counts, i)
          factor_found_counts[i] += 1
        else
          factor_found_counts[i] = 1
        end
      end
    end
  end

  return factor_found_counts
end

# It seems that 5*K*log10(N) samples is often always enough:
N = 10000
K = 3
@time eval_if_lasso_finds_active_factors(N, K, int(5*K*log10(N)), 10)
@time eval_if_lasso_finds_active_factors(N, K, int(5*K*log10(N)), 10; beta_intercept = 0.0)