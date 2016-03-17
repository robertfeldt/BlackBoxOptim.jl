# Based on http://en.wikipedia.org/wiki/Cross-entropy_method

using Distributions

function sample_gaussian(mu, sigma, n)
  distr = Normal(mu, sigma)
  rand(distr, n)
end

# Objective function that we are optimizing
function my_fun1(x)
  exp(-(x-2).^2) + 0.8 * exp(-(x+2).^2)
end

# Cross-entropy using a gaussian proposal distr...
function cross_entropy_method(func;
  max_iterations = 100, num_samples = 100, top_quantile = 0.10,
  std_convergence = 1e-5, mu = 0, sigma = sqrt(100), maximization = true)

  num_top = ceil(Int, num_samples * top_quantile)
  step = 0

  while( step < max_iterations && sigma > std_convergence)
    X = sample_gaussian(mu, sigma, num_samples)
    S = func(X)
    top_X = X[sortperm(S, rev=maximization)[1:num_top]]
    mu = mean(top_X)
    sigma = std(top_X)
    step += 1
  end

  return mu

end
