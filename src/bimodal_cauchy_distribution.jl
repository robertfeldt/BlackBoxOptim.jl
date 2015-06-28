using Distributions

# FIXME implement actual distribution using Distributions.MixtureModel
typealias BimodalCauchy @compat Tuple{Cauchy, Cauchy}

# In the literature Cauchy distributions have been used for sampling the
# f and cr constants used in DE.
function sample_bimodal_cauchy_once(cauchyDistrs::BimodalCauchy, cutoffProb = 0.5)
  index = (rand() < cutoffProb) ? 1 : 2
  rand(cauchyDistrs[index])
end

# When sampling it is common to truncate in either or both ends.
function sample_bimodal_cauchy(cauchyDistrs::BimodalCauchy; cutoffProb = 0.5,
  truncateAbove1 = true, truncateBelow0 = true)
  value = sample_bimodal_cauchy_once(cauchyDistrs, cutoffProb)
  if value > 1.0
    if truncateAbove1
      return 1.0
    else
      return sample_bimodal_cauchy(cauchyDistrs; cutoffProb = cutoffProb, truncateAbove1 = truncateAbove1, truncateBelow0 = truncateBelow0)
    end
  elseif value < 0.0
    if truncateBelow0
      return 0.0
    else
      return sample_bimodal_cauchy(cauchyDistrs; cutoffProb = cutoffProb, truncateAbove1 = truncateAbove1, truncateBelow0 = truncateBelow0)
    end
  else
    return value
  end
end

# For sampling f in DE, bimodal_cauchy(0.65, 0.1, 1.0, 0.1) has been proposed.
# For sampling cr in DE, bimodal_cauchy(0.1, 0.1, 0.95, 0.1) has been proposed.
function bimodal_cauchy(location1, scale1, location2, scale2)
  (Cauchy(location1, scale1), Cauchy(location2, scale2))
end
