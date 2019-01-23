# Compass search as described on page 18 (402) in Kolda2003:
#  Kolda, Tamara G., Robert Michael Lewis, and Virginia Torczon. "Optimization
#  by direct search: New perspectives on some classical and modern methods.",
#  SIAM review 45.3 (2003): 385-482.
function compass_search(f, n; max_fevals = 1e6, delta_tol = 1e-20,
  step_size = 1.0, x = false, known_fmin = false)

  @assert delta_tol > 0
  @assert step_size > delta_tol

  if x == false
    x = randn(n, 1)
  end

  directions = [Matrix{Float64}(I, n,n) -Matrix{Float64}(I, n, n)]

  num_fevals = 0
  fc = fbest = Inf
  candidate = zeros(n, 1)

  while(step_size > delta_tol && num_fevals <= max_fevals)
    # Check all directions to find a better point
    found_better = false
    for direction in 1:(2*n)
      candidate = x + step_size .* directions[:, direction]

      fc = f(candidate)
      num_fevals += 1

      if fc < fbest
        found_better = true
        break
      end

    end

    if found_better
      fbest = fc
      x = candidate

      if fbest == known_fmin
        break
      end
    else
      step_size = step_size / 2
    end
  end

  return x, fbest, num_fevals
end
