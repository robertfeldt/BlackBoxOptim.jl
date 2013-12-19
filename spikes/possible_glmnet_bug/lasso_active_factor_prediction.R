library(glmnet)
N = 10000   # Num factors
K = 3       # Num active factors
active = sample(1:N, K, replace=FALSE)
beta_intercept = 10.0
beta_slope = 1.0
# Lets do positive betas first:
betas = abs(beta_intercept + beta_slope * rnorm(K))
ofunc <- function(x) {
  sum(x[active] * betas)
}
S = round(10*K*log10(N))
xs = matrix(rnorm(N*S), S, N);
y = apply(xs, 1, ofunc)
f = cv.glmnet(xs,y)
pbetas = as.list(coef(f, s="lambda.min"))
# This should return the active factor set +1 since intercept also added.
which( pbetas!=0, arr.ind=TRUE) == ( 1 + sort(c(0, active)) )

# Now lets redo with negative betas:
betas = -betas
ofunc <- function(x) {
  sum(x[active] * betas)
}
y = apply(xs, 1, ofunc)
f = cv.glmnet(xs,y)
pbetas = as.list(coef(f, s="lambda.min"))
# This should return the active factor set +1 since intercept also added.
which( pbetas!=0, arr.ind=TRUE) == ( 1 + sort(c(0, active)) )
