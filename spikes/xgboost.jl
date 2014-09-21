# Testing if we can use XGBoost for surrogate modeling
using XGBoost

N, P, k = (400, 50, 3)
active_coefs = 2 + 10*rand(k) 
coefs = shuffle(vcat(active_coefs, zeros(P-k)))
X = randn(N, P)
Y = X * coefs
testX = randn(N, P)
testY = testX * coefs

mse(y, yhat) = mean((y .- yhat).^2)
made(y, yhat) = mean(abs( 100.0 * (y .- yhat) ./ y ))

num_round = 100
@time bst = xgboost(X, num_round, label = Y, max_depth = 3)
predY = predict(bst, testX);
mse(testY, predY)
made(testY, predY)
