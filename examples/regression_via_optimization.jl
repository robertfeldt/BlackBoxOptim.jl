# Lets use black-box optimization to do regression!
#
# The basic idea is that a regression problem is a problem of deciding
# on which set of coefficients (beta) that minimizes an objective function
# that compares the dependent variables y to their independent variables x.
# When implementing this x will be a matrix d*n of floats where d is the dimension
# of the problem and n is the number of examples/cases. Thus y is an array/vector
# 1*n of floats and the coefficients is an array/vector (d+1)*1 of coefficients
# where the first is the intercept (beta0) while the rest are the coefficients.
#
# We get different types of regression depending on which objective function
# we select. Common to many of them is that they first calculate the deviances
# between the model and the actual values the model should predict:
function discrepancies(beta, x, y)
  beta0 = beta[1]
  dims = length(beta) - 1
  y - ( (ones(1, size(x, 2)) * beta0) .+ sum(beta[2:end] .* x, 1) )
end

# Given this setup we can now create an objective function for 
# Ordinare Least Squares (OLS) regression. This is actually just the L2 norm:
function ols_regression_objective(beta, x, y)
  norm(discrepancies(beta, x, y), 2)
end

# And the L1 norm gives Least Absolute Deviations (LAD) (aka Robust regression) regression:
function lad_regression_objective(beta, x, y)
  norm(discrepancies(beta, x, y), 1)
end

# And we can do LASSO and Ridge regression by adding a penalty on large coefficients, but 
# note that these also take a lambda constant:
function regularized_regression_objective(lambda, beta, x, y, p = 2, q = 1)
  norm(discrepancies(beta, x, y), p) + lambda * norm(beta[2:end], q)
end

# By selecting p = 2 and q = 1 we get LASSO regression:
function lasso_regression_objective(lambda, beta, x, y)
  regularized_regression_objective(lambda, beta, x, y, 2, 1)
end

# and by selecting p = 2 and q = 2 we get Ridge regression:
function ridge_regression_objective(lambda, beta, x, y)
  regularized_regression_objective(lambda, beta, x, y, 2, 2)
end

# Ok, we're ready to do some regression. We need some data, lets start with
# a way to generate data from models specified as julia functions. Lets start
# with a simple model of 3 vars: X1 + 2*X2 - X3
function m1(x)
  x[1,:] + 2*x[2,:] - x[3,:]
end

# and we generate random input for it:
x1 = rand(3, 100)
# and calc the outputs:
y1 = m1(x1)

# We can now search for the coefficients that minimize the OLS objective
# with a black-box optimization search like so (we allow coefficients to have
# a min value of -5.0 and a max value of 5.0, and the search is for 4 
# coefficients, one intercept and three for each of the values of x):
using BlackBoxOptim
ols_bestfit, ols_error = bboptimize((betas) -> ols_regression_objective(betas', x1, y1), (-10.0, 10.0); dimensions = 4, iterations = 2e4)

# But the really nice thing is that we can easily consider other objectives such as the LAD:
lad_bestfit, lad_error = bboptimize((betas) -> lad_regression_objective(betas', x1, y1), (-10.0, 10.0); dimensions = 4, iterations = 2e4)

# For regularized regression we can optimize for different values of lambda so
# create a wrapper function that handles this:
function regularized_opt(lambda, func)
  bboptimize((betas) -> func(lambda, betas', x1, y1), 
    (-10.0, 10.0); dimensions = 4, iterations = 2e4)
end

lasso_bestfit1, lasso_error1 = regularized_opt(1, lasso_regression_objective)
lasso_bestfit2, lasso_error2 = regularized_opt(2, lasso_regression_objective)
lasso_bestfit3, lasso_error3 = regularized_opt(3, lasso_regression_objective)

ridge_bestfit1, ridge_error1 = regularized_opt(1, ridge_regression_objective)
ridge_bestfit2, ridge_error2 = regularized_opt(2, ridge_regression_objective)
ridge_bestfit3, ridge_error3 = regularized_opt(3, ridge_regression_objective)

# Conclusion: With black-box optimization you can easily fit regression models
# from very different paradigms without having to implement very much code.