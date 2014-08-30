include("../src/frequency_adaptation.jl")

# Combine Loschilov's ACD with Glasmachers ACF-CD. Initial testing indicated it is 3-20 times
# worse than just ACD. I guess the ACD needs new fitness evals in each direction of the transformation
# matrix in order to reliably update it. So maybe only use ACF when increasing the population size, i.e.
# ensure that all directions have been tried at least once and then add additional once from the frequency 
# selected directions?
function acf_adaptive_coordinate_descent(fitnessfct, P, Xmin, Xmax;
  MaxEval = 1e4 * P, 
  stopfitness = 1e-8, 
  howOftenUpdateRotation = 1,
  numfevals = 0,
  kSuccess = 2.0,
  kUnsuccess = 0.5,
  alwaysSampleTwoPoints = false, # sample 2 points per coordinate regardless if 1st is improvement
  nonImprovementBudget = 1       # Percent of fevals we allow before restarting if not enough improvement
  )

  # Initialize mean and sigma
  xmins = Xmin * ones(P)
  xmaxs = Xmax * ones(P)
  xmean = xmins .+ rand(P) * (Xmax - Xmin)
  sigma = (xmaxs .- xmins) / 4 # Step size per principal component (NOT per dimension!)

  bestFit = 1e+100 # Start from a very high fitness to ensure we will find lower... Fix smarter...

  mu = P # Use the P best candidates among the 2*P in each update of the adaptive encoding...
  ae = AdaptiveEncoding(mu, P)
  # To compare to CD without adaptive encoding instead use:
  #ae = NoEncoding(mu, P)
  #alwaysSampleTwoPoints = true

  iter = 0
  somebetter = false

  # There can be 2*P or 2*P+1 candidate points tested before we update the encoding,
  # since we only add the 2nd point if the first was not better.
  allx = zeros(P, 2*P+1) # Each column is a candidate and we will save 2 candidates per component/iteration in the population
  allf = zeros(2*P+1)    # One fitness value per candidate
  currentpos = 0

  # Restart the search from a new random position unless the fitness has improved more than 2*stopfitness
  # in last nonImprovementBudget % of fevals.
  lastconvergencecheck = numfevals
  lastcheckbestfit = bestFit

  # We will adaptively select among the P directions
  freqadapter = FrequencyAdapter(P)

  while ((numfevals < MaxEval) && (bestFit > stopfitness))
    ix = next(freqadapter) # Get the next direction according to the adapted frequencies

    # Sample two candidate solutions and calc their fitness
    dx = sigma[ix] * ae.B[:,ix] # shift along ix'th principal component, the computational complexity is linear

    # Randomly take the minus or plus direction first.
    direction = (rand() < 0.5) ? 1 : (-1)

    x1 = xmean + direction * dx             # first point to test along ix'th principal component
    Fit1 = fitnessfct(x1)
    update!(freqadapter, ix, (bestFit - Fit1) / bestFit)
    numfevals += 1

    # Update archive with new point
    # posx1 = ix*2 - 1
    currentpos += 1
    allx[:, currentpos] = x1
    allf[currentpos] = Fit1

    # Found better point?
    foundbetter = false
    if (Fit1 < bestFit)
      bestFit = Fit1
      fevalsatbest = numfevals
      xmean = x1
      foundbetter = true
    end

    # Only do second point if first point was not better (unless parameter states we should always smaple two).
    if !foundbetter || alwaysSampleTwoPoints
      x2 = xmean - direction * dx             # second point to test is symmetric to the first one on the same ix'principal component
      Fit2 = fitnessfct(x2)
      update!(freqadapter, ix, (bestFit - Fit2) / bestFit)
      numfevals += 1
      # Update archive with new point
      currentpos += 1
      # posx2 = ix*2
      allx[:, currentpos] = x2
      allf[currentpos] = Fit2

      if (Fit2 < bestFit)
        bestFit = Fit2
        fevalsatbest = numfevals
        xmean = x2
        foundbetter = true
      end
    end

    # Adapt step-size sigma depending on the success/unsuccess of the last search step.
    # Increase step size if found a better candidate, decrease if not.
    if foundbetter
      somebetter = true
      sigma[ix] *= kSuccess
    else
      sigma[ix] *= kUnsuccess
    end

    # Update the encoding every time we have tested at least 2*P candidate points and at least one
    # better candidate found.
    if currentpos >= 2*P
      if somebetter
        somebetter = false
        fitnessperm = sortperm(allf[1:currentpos])
        population = allx[:,fitnessperm]
        update!(ae, population, howOftenUpdateRotation)
      end

      # Check for convergence/nonimprovement every now and then and restart if so.
      if (numfevals - lastconvergencecheck) >= (MaxEval * nonImprovementBudget / 100.0)
        # We make a simple prediction of the final fitness given the slope between last
        # check and this. Restart if the predicted is too far from stopfitness.
        pfitness = predicted_final_fitness(numfevals, lastconvergencecheck, MaxEval, bestFit, lastcheckbestfit)
        if pfitness > 1000*stopfitness
          println("Restart!")
          # Recursive call to restart but express how many function evals we have already done...
          xm2, bf2, nf2 = acf_adaptive_coordinate_descent(fitnessfct, P, Xmin, Xmax;
            MaxEval = MaxEval, stopfitness = stopfitness, 
            howOftenUpdateRotation = howOftenUpdateRotation, 
            numfevals = numfevals)
          if bf2 < bestFit
            return (xm2, bf2, nf2)
          else
            return (xmean, bestFit, nf2)
          end
        end
        lastcheckbestfit = bestFit
        lastconvergencecheck = numfevals
      end

      currentpos = 0
    end

    iter += 1
    if rand() < 0.01
      pf = predicted_final_fitness(numfevals, lastconvergencecheck, MaxEval, bestFit, lastcheckbestfit)
      println("$(numfevals): bestfit = $(bestFit) (predicted: $(pf))")
    end
  end

  return (xmean, bestFit, numfevals)
end

function predicted_final_fitness(currentfevals, lastfevals, maxfevals, currentfitness, lastfitness)
  slope = (currentfitness - lastfitness) / (currentfevals - lastfevals)
  currentfitness + slope * (maxfevals - currentfevals)
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

# Use this when no adaptive encoding is used, i.e. the ACD algorithm is a normal CD.
type NoEncoding
  B
  NoEncoding(mu, P) = begin
    new(eye(P, P))
  end
end

function update!(e::NoEncoding, population::Matrix{Float64}, howOftenUpdateRotation = 1)
  # Do nothing since we never update anything...
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
    pop = population[:,1:ae.mu]
  else
    pop = population
  end
  ae.xmean = pop * ae.weights  # Line 9

  # Line 10: Calculate z0 but check for the denominator being zero...
  xdiff = ae.xmean .- xold
  denom = sum((ae.invB * xdiff).^2)
  if denom == 0
    z0 = zeros(ae.P)
  else
    z0 = (sqrt(ae.P) / sqrt(denom)) * xdiff
  end
  
  # Line 11&13: Calc zi's and then Cmu.
  #ae.cmu = 0.0
  if ae.cmu != 0.0
    xis_minus_xold = broadcast(-, pop, xold)
    denoms = sum((ae.invB * xis_minus_xold).^2, 1)
    zis = broadcast(*, (sqrt(ae.P) ./ sqrt(denoms)), xis_minus_xold)
    Cmu = broadcast(*, ae.weights', zis * zis')
  else
    Cmu = 1 # dummy just so it is defined
  end

  # Line 12: Update p 
  ae.pc = (1-ae.cp) * ae.pc .+ sqrt(ae.cp*(2-ae.cp)) * z0

  # Line 14: Update C.
  ae.C = (1-ae.c1 - ae.cmu) * ae.C .+ ae.c1 * (ae.pc * ae.pc') .+ ae.cmu * Cmu

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
