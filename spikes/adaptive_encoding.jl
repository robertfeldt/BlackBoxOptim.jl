abstract type CoordinateTransform end

# Based on the pseudo-code on page 156 in Loschchilov's PhD Thesis from 2013.
mutable struct AdaptiveEncoding <: CoordinateTransform
  n::Int
  sqrtn::Float64
  initialized::Bool
  c_p::Float64
  c_mu::Float64
  c_1::Float64
  p::Array{Float64, 2}
  C::Array{Float64, 2}
  B::Array{Float64, 2}
  m::Array{Float64, 2}
  m_prev::Array{Float64, 2}

  AdaptiveEncoding(mu, n) = begin
    c_p = 1.0 / sqrt(n)
    c_mu = c_1 = 0.5 / n
    p = zeros(n, 1)
    C = Matrix{Float64}(I, n, n)
    B = Matrix{Float64}(I, n, n)
    new(n, sqrt(n), false, c_p, c_mu, c_1,
      p, C, B, zeros(n, 1), zeros(n, 1))
  end
end

# Update the coordinate transform given a new set of points (x) and their
# utilities/fitnesses (us). Based on
# http://en.wikipedia.org/wiki/CMA-ES#Example_code_in_MATLAB.2FOctave
function update_transform!(ae::AdaptiveEncoding, x::Array{Float64, 2}, us::Array{Float64,2})
  # Init m on the first call.
  if !ae.initialized
    ae.m = x * us
    ae.initialized = true
    # The other variables were already updated in the constructor...
    return ae.B
  end

  ae.m_prev = ae.m
  ae.m = x * us
  Binv = inv(ae.B)

  diff_m = ae.m - ae.m_prev
  z0 = (ae.sqrtn / norm(Binv * diff_m)) * diff_m

  # Update evolution path.
  ae.p = (1.0 - ae.c_p) * ae.p + sqrt(ae.c_p * (2 - ae.c_p)) * z0

  # Skip the rank-mu update when you need to save time.
  # Only do this randomly, at times, then, to save time?
  xis_minus_m_prev = broadcast(-, xi, m_prev)
  zis = (ae.sqrtn / norm(Binv * xis_minus_m_prev)) * xis_minus_m_prev
  C_mu = zeros(ae.n, ae.n)
  # Should be some simpler way to achieve this... Maybe C_mu = zis * diag(us) * zis'
  for i in 1:ae.n
    z = zis[:,i]
    C_mu += (us[i] * (z * z'))
  end

  # Now update the covar matrix. Skip the c_mu parts if not doing the rank-u update...
  ae.C = (1 - ae.c_1 - ae.c_mu) * ae.C + (ae.c_1 * ae.p * ae.p') + ae.c_mu * C_mu

  # Should we enforce symmetry in C?
  # tc = triu(ae.C)
  # ae.C = tc + tc' - Matrix{Float64}(I, ae.C, ae.C) * diag(ae.C)

  # Now B is the square root of the covar matrix, but faster to do eig decomposition.
  eigvalues, eigvectors = eig(ae.C)
  diagD = diag(eigvectors)
  ae.B = eigvalues .* repmat(sqrt(diagD)', ae.n, 1)

  return ae.B
end
