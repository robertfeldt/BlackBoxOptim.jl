function adaptive_coordinate_descent(fitnessfct, P, Xmin, Xmax;
  MaxEval = 1e5 * P, 
  stopfitness = 1e-10, 
  howOftenUpdateRotation = 1,
  numfevals = 0,
  kSuccess = 2.0,
  kUnsuccess = 0.5,
  c1 = 0.5 / P,
  cmu = 0.5 / P,
  nonImprovementBudget = 1 # Percent of remaining budget that can be used without fitness improvement before restart
  )

  @assert howOftenUpdateRotation >= 1

  # Initialize mean and sigma
  xmins = Xmin * ones(P)
  xmaxs = Xmax * ones(P)
  xmean = xmins .+ rand(P) * (Xmax - Xmin)
  sigma = (xmaxs .- xmins) / 4 # Step size per principal component (NOT per dimension!)

  bestFit = 1e+100 # Start from a very high fitness to ensure we will find lower... Fix smarter...

  mu = P # Use the P best candidates among the 2*P in each update of the adaptive encoding...
  ae = AdaptiveEncoding(mu, P)

  iter = 0
  xperm = shuffle(collect(1:P)) # Use them in random order...
  qix = 0
  somebetter = false
  allx = zeros(P, 2*P) # Each column is a candidate and we will save 2 candidates per component/iteration in the population
  allf = zeros(2*P)    # One fitness value per candidate

  while ((numfevals < MaxEval) && (bestFit > stopfitness))
    qix = 1 + iter % P
    ix = xperm[qix]        # will perform the search along this ix'th principal component

    # Sample two candidate solutions and calc their fitness
    dx = sigma[ix] * ae.B[:,ix] # shift along ix'th principal component, the computational complexity is linear
    x1 = xmean - dx             # first point to test along ix'th principal component
    x2 = xmean + dx             # second point to test is symmetric to the first one on the same ix'principal component
    Fit1 = fitnessfct(x1)
    Fit2 = fitnessfct(x2)
    numfevals += 2

    # Found better point?
    foundbetter = false
    if (Fit1 < bestFit)
      bestFit = Fit1
      fevalsatbest = numfevals
      xmean = x1
      foundbetter = true
    end
    if (Fit2 < bestFit)
      bestFit = Fit2
      fevalsatbest = numfevals
      xmean = x2
      foundbetter = true
    end

    # Adapt step-size sigma depending on the success/unsuccess of the last search step.
    # Increase step size if found a better candidate, decrease if not.
    if foundbetter
      somebetter = true
      sigma[ix] *= kSuccess
    else
      sigma[ix] *= kUnsuccess
    end

    # Update archive 
    posx1 = ix*2 - 1
    posx2 = ix*2
    allx[:, posx1] = x1
    allf[posx1] = Fit1
    allx[:, posx2] = x2
    allf[posx2] = Fit2

    # Update the encoding every P steps if at least one better candidate found among last P
    # steps.
    if qix == P
      xperm = shuffle(collect(1:P))
      if somebetter
        somebetter = false
        fitnessperm = sortperm(allf)
        population = allx[:,fitnessperm]
        update!(ae, population, howOftenUpdateRotation)
      else # A full iteration without any improvement; check if we need to restart
        # If no improvement in a long time we restart
        if (numfevals - fevalsatbest) > max(20*P, (MaxEval - fevalsatbest) * nonImprovementBudget / 100.0)
          println("Restart!")
          # Recursive call to restart but express how many function evals we have already done...
          xm2, bf2, nf2 = adaptive_coordinate_descent(fitnessfct, P, Xmin, Xmax;
            MaxEval = MaxEval, stopfitness = stopfitness, 
            howOftenUpdateRotation = howOftenUpdateRotation, 
            numfevals = numfevals)
          if bf2 < bestFit
            return (xm2, bf2, nf2)
          else
            return (xmean, bestFit, nf2)
          end
        end
      end
    end

    iter += 1
    if rand() < 0.02
      println("$(numfevals): bestfit = $(bestFit) ($(numfevals - fevalsatbest), $(max(20*P, (MaxEval - fevalsatbest) * nonImprovementBudget / 100.0)))")
    end
  end

  return (xmean, bestFit, numfevals)
end

function generate_random_orthogonal_transformation(P)
  B = eye(P,P)
  for i = 1:P
    v = randn(P,1)
    while norm(v) < 1e-2
      v = randn(P,1)
    end
    for j = 1:i-1
      v = v - (v'*B[:,j]) .* B[:,j]
    end
    B[:,i] = v / norm(v)
  end
  return B
end

type AdaptiveEncoding
  mu
  P
  B
  Bo
  invB
  weights
  mucov
  alpha_p
  c1
  cmu
  cp
  pc
  pcmu
  xmean
  C
  Cold
  diagD
  ps
  iter
  AdaptiveEncoding(mu::Int64, P::Int64, B::Matrix{Float64}, population = randn(P, mu)) = begin
    # We assume B is orthogonal!
    weights = ones(mu) / mu # Can also try non-uniform weights here...
    new(mu, P, B, B, B',
      weights,
      mu, # mucov
      1, # alpha_p
      0.5 / P, # c1
      0.5 / P, # cmu
      1 / sqrt(P), # cp
      zeros(P), # pc
      zeros(P), # pcmu
      population * weights, # xmean
      eye(P), # C
      eye(P), # Cold
      ones(P), # diagD
      0, # ps
      0  # iter
    )
  end
  AdaptiveEncoding(mu::Int64, P::Int64) = begin
    AdaptiveEncoding(mu, P, generate_random_orthogonal_transformation(P))
  end
end

function update!(ae::AdaptiveEncoding, population::Matrix{Float64}, howOftenUpdateRotation = 1)
  ae.iter += 1

  #
  # Update the adaptive ancoding, see lines 8-14 on Algorithm 2 of the Loschilov ACD paper.
  #

  xold = ae.xmean # Line 8

  # Line 9 but ensure only first mu candidates of population is used => they should be the best ones
  # in fitness-sorted order.
  if size(population, 2) != ae.mu
    ae.xmean = population[:,1:ae.mu] * ae.weights  # Line 9
  else
    ae.xmean = population * ae.weights
  end

  # Line 10: Calculate z0 but check for the denominator being zero...
  xdiff = ae.xmean .- xold
  denom = sum((ae.invB * xdiff).^2)
  if denom == 0
    z0 = zeros(ae.P)
  else
    z0 = (sqrt(ae.P) / sqrt(denom)) * xdiff
  end
  
  # Line 11&13: Calc zi's and then Cmu. TBD.

  # Line 12: Update p 
  ae.pc = (1-ae.cp) * ae.pc .+ sqrt(ae.cp*(2-ae.cp)) * z0

  # Line 14: Update C. Cmu part not yet implemented.
  ae.C = (1-ae.c1) * ae.C + ae.c1 * (ae.pc * ae.pc')

  if ae.iter % howOftenUpdateRotation == 0 || ae.iter <= 2
    # Go back to old covar matrix if latest has NaN or Inf element
    if any(isnan(ae.C)) || any(isinf(ae.C))
      ae.C = ae.Cold
    end

    # Ensure covar matrix is symmetric by copying the upper triangular part to the lower tri part.
    ae.C = (triu(ae.C) + triu(ae.C, 1)')

    # Now do the eigen decomposition. This is the costly operation...
    eigenvalues, ae.Bo = eig(ae.C)

    # Limit the condition value of C and control for extreme eigen values.
    cond = 1e12
    if minimum(eigenvalues) <= 0
      eigenvalues[eigenvalues .< 0] = 0
      tmp = maximum(eigenvalues) / cond
      ae.C = ae.C .+ tmp * eye(ae.P)
      eigenvalues += tmp * ones(ae.P,1)
    end
    if maximum(eigenvalues) > cond*minimum(eigenvalues)
      tmp = maximum(eigenvalues)/cond - minimum(eigenvalues)
      ae.C = ae.C .+ tmp * eye(ae.P)
      eigenvalues += tmp*ones(ae.P,1)
    end
    if (minimum(eigenvalues) <= 0)
      # Should never happen since we take care of this in the limiting code above...
      throw("Negative eigenvalues => skipping update")
    end

    try
      ae.diagD = sqrt(eigenvalues) 
      # Use Cholesky instead??
      ae.B = ae.Bo * diagm(ae.diagD) # This is sometimes a matrix rather than a vector. Strange!
      ae.invB = diagm(1 ./ ae.diagD) * ae.Bo'
      ae.Cold = ae.C
    catch
      # Reset to old covar matrix with a random perturbation to try to get out of the problems
      ae.C = ae.Cold .+ 0.01 * rand(ae.P, ae.P)
    end
  end
end
