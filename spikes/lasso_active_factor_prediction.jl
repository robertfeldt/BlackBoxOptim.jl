# Idea: Use LASSO for screening the active factors in an optimization
# problem.

using GLMNet

# Say that we want to optimize a function which has many factors but where
# only very few of them are active in calculating the outcome.
N = 1000
K = 3
active_factors = shuffle(collect(1:N))[1:K]
betas = 5.0 + 5.0 * randn(K, 1)
function ofunc(x)
  active = x[active_factors]
  sum(active .* betas)
end

# No lets sample
S = int(N/20)
xs = randn(N, S)

# And calc the function values for them:
y = [ofunc(xs[:,i]) for i in 1:S]

# Now do lasso regression:
# path = glmnet(xs', y)

function active_factors_from_lasso(xs, y)
  # Find the best lambda using cross-validation:
  cv = glmnetcv(xs, y)
  num_steps = size(cv.path.betas, 2)
  res = Any[]
  for(i in 1:num_steps)
    pbetas = cv.path.betas[:, i]
    pbwi = [(pbetas[i], i) for i in 1:length(pbetas)]
    push!(res, filter((p) -> p[1] > 0.01, pbwi))
  end
  return res
end

lasso_factors = active_factors_from_lasso(xs', y)
active_factors