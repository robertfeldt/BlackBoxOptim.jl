# Lets use black-box optimization to do regression!
#
# To run this, just do:
#   julia regression_via_optimization.jl
# from a command line after having installed BlackBoxOptim.
#
# The basic idea is that a regression problem is a problem of deciding
# which set of coefficients (betas) minimize an objective function
# that compares the dependent variable y to a model based on the betas and the
# independent variables x.
#
# When implementing this, x will be a matrix d*n of floats where d is the dimension
# of the problem and n is the number of examples/cases. Thus y is an array/vector
# 1*n of floats and the coefficients is an array/vector (d+1)*1 of coefficients
# where the first is the intercept (beta0) while the rest are the coefficients.
#
# We get different types of regression depending on which objective function
# we select. Common to many of them is that they first calculate the deviances
# between the model and the actual values the model should predict:
function discrepancies(beta, x, y)
  beta0 = beta[1]
  ncols = size(x, 2)
  y - ( (ones(1, ncols) * beta0) .+ sum(beta[2:end] .* x, 1) )
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
ols_bestfit, ols_error = bboptimize((betas) -> ols_regression_objective(betas', x1, y1),
  (-10.0, 10.0); dimensions = 4, iterations = 2e4)

# But the really nice thing is that we can easily consider other objectives such as the LAD:
lad_bestfit, lad_error = bboptimize((betas) -> lad_regression_objective(betas', x1, y1),
  (-10.0, 10.0); dimensions = 4, iterations = 2e4)

# For regularized regression we can optimize for different values of lambda so
# create a wrapper function that handles this:
function regularized_opt(lambda, func, x, y, dims, its = 2e4)
  bboptimize((betas) -> func(lambda, betas', x, y),
    (-10.0, 10.0); dimensions = dims, iterations = its)
end

lasso_bestfit1, lasso_error1 = regularized_opt(1, lasso_regression_objective, x1, y1, 4)
lasso_bestfit2, lasso_error2 = regularized_opt(2, lasso_regression_objective, x1, y1, 4)
lasso_bestfit3, lasso_error3 = regularized_opt(3, lasso_regression_objective, x1, y1, 4)

ridge_bestfit1, ridge_error1 = regularized_opt(1, ridge_regression_objective, x1, y1, 4)
ridge_bestfit2, ridge_error2 = regularized_opt(2, ridge_regression_objective, x1, y1, 4)
ridge_bestfit3, ridge_error3 = regularized_opt(3, ridge_regression_objective, x1, y1, 4)

# Now lets create some support functions for printing models nicely.
linear_terms(num) = [@sprintf(" * X%d", i) for i in 1:num]
squared_terms(num) = [@sprintf(" * X%d^2", i) for i in 1:num]
linsq_terms(num) = vcat(linear_terms(num), squared_terms(num))

function sprint_predicted_model(bestfit, terms = nothing, skipIfLower = 1e-5)
  if terms == nothing
    terms = linear_terms(length(bestfit)-1)
  end
  if length(terms) < length(bestfit)
    terms = vcat([""], terms) # Put an empty term first which corresponds to intercept
  end

  signstr(value) = (value < 0.0) ? (@sprintf(" - %.3f", abs(value))) : (@sprintf(" + %.3f", value))

  elems = Any[]
  first_push = true
  for i in 1:length(bestfit)
    if abs(bestfit[i]) > skipIfLower
      str = join([signstr(bestfit[i]), terms[i]])
      if first_push
        # Strip away leading sign since this it the first push
        push!(elems, str[3:end])
        first_push = false
      else
        push!(elems, str)
      end
    end
  end

  join(elems)
end

# Let's try a model which involves a squared terms:
#   X1 + 4.13*X2*X2 - 3.14*X3
function m2(x)
  x[1,:] + (4.13 * x[2,:].^2) - (3.14 * x[3,:])
end

# and we generate random input for it:
x2 = rand(3, 100)
# and calc the outputs:
y2 = m2(x2)

# Before we regress we need to encode our beliefs about the general
# structure of the model. Let's say we believe there are squared terms but we
# do not know which ones. So we add one squared term per independent variable:
x2m = zeros(3+3, 100)
x2m[1:3,:] = x2
x2m[4,:] = x2[1,:].^2
x2m[5,:] = x2[2,:].^2
x2m[6,:] = x2[3,:].^2

# With this we can fit models:
m2_ols_bestfit, m2_ols_error = bboptimize((betas) -> ols_regression_objective(betas', x2m, y2),
  (-10.0, 10.0); dimensions = 1+3+3, iterations = 5e4)
m2_lad_bestfit, m2_lad_error = bboptimize((betas) -> lad_regression_objective(betas', x2m, y2),
  (-10.0, 10.0); dimensions = 1+3+3, iterations = 5e4)
m2_lasso_bestfit1, m2_lasso_error1 = regularized_opt(1, lasso_regression_objective, x2m, y2, 7, 5e4)
m2_lasso_bestfit2, m2_lasso_error2 = regularized_opt(2, lasso_regression_objective, x2m, y2, 7, 5e4)
m2_lasso_bestfit3, m2_lasso_error3 = regularized_opt(3, lasso_regression_objective, x2m, y2, 7, 5e4)
m2_ridge_bestfit1, m2_ridge_error1 = regularized_opt(1, ridge_regression_objective, x2m, y2, 7, 5e4)
m2_ridge_bestfit2, m2_ridge_error2 = regularized_opt(2, ridge_regression_objective, x2m, y2, 7, 5e4)
m2_ridge_bestfit3, m2_ridge_error3 = regularized_opt(3, ridge_regression_objective, x2m, y2, 7, 5e4)

# And now lets print our models nicely so user can see the results...

println("Model1 = 1.000 * X1 + 2.000 * X2 - 1.000 * X3")
println("OLS best fit: ", sprint_predicted_model(ols_bestfit))
println("LAD best fit: ", sprint_predicted_model(lad_bestfit))
println("LASSO best fit, lambda = 1: ", sprint_predicted_model(lasso_bestfit1))
println("LASSO best fit, lambda = 2: ", sprint_predicted_model(lasso_bestfit2))
println("LASSO best fit, lambda = 3: ", sprint_predicted_model(lasso_bestfit3))
println("Ridge best fit, lambda = 1: ", sprint_predicted_model(ridge_bestfit1))
println("Ridge best fit, lambda = 2: ", sprint_predicted_model(ridge_bestfit2))
println("Ridge best fit, lambda = 3: ", sprint_predicted_model(ridge_bestfit3))

terms = linsq_terms(3)
println("")
println("Model2 = 1.000 * X1 - 3.140 * X3 + 4.130 * X2^2 ")
println("OLS best fit: ", sprint_predicted_model(m2_ols_bestfit, terms))
println("LAD best fit: ", sprint_predicted_model(m2_lad_bestfit, terms))
println("LASSO best fit, lambda = 1: ", sprint_predicted_model(m2_lasso_bestfit1, terms))
println("LASSO best fit, lambda = 2: ", sprint_predicted_model(m2_lasso_bestfit2, terms))
println("LASSO best fit, lambda = 3: ", sprint_predicted_model(m2_lasso_bestfit3, terms))
println("Ridge best fit, lambda = 1: ", sprint_predicted_model(m2_ridge_bestfit1, terms))
println("Ridge best fit, lambda = 2: ", sprint_predicted_model(m2_ridge_bestfit2, terms))
println("Ridge best fit, lambda = 3: ", sprint_predicted_model(m2_ridge_bestfit3, terms))

# Conclusion: With black-box optimization you can easily fit regression models
# from very different paradigms without having to implement very much code.
