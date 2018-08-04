using GLMNet
N = 10000
K = 3
active = shuffle(collect(1:N))[1:K]
beta_intercept = 10.0
beta_slope = 1.0
betas = abs(beta_intercept + beta_slope * randn(K, 1))
ofunc(x) = sum(x[active] .* betas)
S = int(10*K*log10(N));
xs = randn(S, N);
y = [ofunc(xs[i,:]) for i in 1:S];
cv = glmnetcv(xs, y; standardize = false)
pbetas = cv.path.betas[:, argmin(cv.meanloss)]
pbwi = [(pbetas[i], i) for i in 1:length(pbetas)]
filter((p) -> abs(p[1]) > 0.01, pbwi)

betas = -betas
ofunc(x) = sum(x[active] .* betas)
y = [ofunc(xs[i,:]) for i in 1:S];
cv = glmnetcv(xs, y; standardize = false)
pbetas = cv.path.betas[:, argmin(cv.meanloss)]
pbwi = [(pbetas[i], i) for i in 1:length(pbetas)]
filter((p) -> abs(p[1]) > 0.01, pbwi)
