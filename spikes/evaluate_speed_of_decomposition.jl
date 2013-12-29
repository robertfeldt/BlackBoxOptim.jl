# Speed test for decomposing matrices in Julia.

function speed_test(func, orderN = 1000)
  sizes = Any[]
  for i in 0:int(log10(orderN))
    push!(sizes, 1 * 10^i)
    push!(sizes, 2 * 10^i)
    push!(sizes, 5 * 10^i)
  end

  y = zeros(length(sizes))
  for i in 1:length(sizes)
    n = sizes[i]
    A = randn(n, n)
    A = A'*A
    print("Testing with n = $(n)")
    tic()
    res = func(A)
    y[i] = toq()
    println(", time = $(y[i]), ")
  end

  return sizes, y
end

n_eig, t_eig = speed_test(eig, 1000)
n_chol, t_chol = speed_test(chol, 1000)

# Linear model
lin_eig = linreg(t_eig, n_eig)