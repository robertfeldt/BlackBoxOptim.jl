# Basic experimentation with SparseCholeskyCovarSampler indicated
# that it is hard to not run into exceptions such as non-definite matrices etc
# when the dimensionality increases. Since high dimensionality is the motivation
# for using a sparse variant I have not experimented further with this...

# This is not a general solution but it works for specifically my use case and
# when the operator is *
function Base.broadcast{Tv,Ti}(op, v::Array{Tv,2}, A::SparseMatrixCSC{Tv,Ti})
  I, J = findn(A)
  V = zeros(nnz(A))
  vn, vm = size(v)
  if vn >= 1 && vm == 1
    for(l in 1:nnz(A))
      row = I[l]
      V[l] = op(v[row], A.nzval[l])
    end
  elseif vn == 1 && vm >= 1
    for(l in 1:nnz(A))
      col = J[l]
      V[l] = op(v[col], A.nzval[l])
    end
  else
    throw(ArgumentError("invalid dimensions"))
  end
  sparse(I, J, V)
end

type SparseCholeskyCovarSampler <: CovarianceMatrixSampler
  C::SparseMatrixCSC{Float64,Int64}
  sqrtC::SparseMatrixCSC{Float64,Int64}

  SparseCholeskyCovarSampler(n) = begin
    new(speye(n,n), speye(n,n))
  end
end

function update_covariance_matrix!(cms::SparseCholeskyCovarSampler, delta, a)
  C = a * cms.C + (1 - a) * delta
  cms.C = C # triu(C) + triu(C,1)' # Ensure C is symmetric. Should not be needed, investigate...
end

function decompose!(cms::SparseCholeskyCovarSampler)
  try
    #println("droptol!(C)")
    #cms.C = Base.droptol!(cms.C, 1e-4)
    #println("nnz(C) = $(nnz(cms.C)), min = $(minimum(cms.C.nzval)), max = $(maximum(cms.C.nzval))")
    #println("cholfact(C)")
    t = cholfact(cms.C)
    #println("sparse(t)")
    s = sparse(t)
    #println("s'")
    sqrtC = s'
    #println("droptol!(sqrtC)")
    #cms.sqrtC = Base.droptol!(sqrtC, 1e-4)
    #show(cms.sqrtC)
  catch error
    # We don't update if there is some problem
    println("ERROR: ", error)
  end
end

function multivariate_normal_sample(cms::SparseCholeskyCovarSampler, n, m)
  #if n*m > 100
    # We select a random density in [0.10, 0.60]
  #  density = 0.10 + 0.60 * rand()
  #else
  #  density = 0.90
  #end
  cms.sqrtC * randn(n, m)
end
