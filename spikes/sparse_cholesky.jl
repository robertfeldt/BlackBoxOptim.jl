# In covariance adaptation we want to use cholesky factorization instead
# of a full eigendecomposition. Can we speed this up with sparse matrices in
# Julia?

# Objective function to be maximized.
function my_fun1(x)
  exp(-(x-2).^2) + 0.8 * exp(-(x+2).^2)
end

Ns = [10, 100, 1000, 2000]
num_reps = 10
Ds = [0.01, 0.10, 0.50]

function non_sparse_sqrtm(C, randsamples)
  sqrtC = sqrt(C)
  s = sqrtC * randsamples
  C + (s * s')
end

function sparse_cholfact(C, randsamples)
  sqrtC = sparse(cholfact(C))
  s = sqrtC' * randsamples
  C = C + (s * s')
end

for n in Ns
  # Non-sparse, sqrtm
  C = Matrix{Float64}(I, n, n)
  if n <= 2000
    tic()
    for rep in 1:num_reps
      C = non_sparse_sqrtm(C, randn(n, n))
    end
    t = toq()
    tpdn = t/(num_reps*n)
    println("Non-sparse, sqrtm, n = $(n), time = $(t/num_reps), time = $(tpdn)")
  end

  for d in Ds
    # Sparse, cholfact, randn
    C = speye(n)
    tic()
    for rep in 1:num_reps
      C = sparse_cholfact(C, sprandn(n, n, d))
    end
    t = toq()
    tpds = t/(num_reps*n)
    println("Sparse, cholfact, d = $(d), n = $(n), time = $(t/num_reps), time = $(tpds), $(tpdn/tpds) times faster")
  end
end
