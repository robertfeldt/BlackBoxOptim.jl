module Utils

# Latin hypercube sampling of values from given mins and maxs, per column.
function latin_hypercube_sampling(mins, maxs, numSamples)
  dims = length(mins)
  result = zeros(numSamples, dims)
  for(i in 1:dims)
    interval_len = (maxs[i] - mins[i]) / numSamples
    interval_starts = linspace(mins[i], maxs[i], numSamples+1)[1:numSamples]
    samples_for_dim = interval_starts + interval_len * rand(numSamples, 1)
    result[:,i] = samples_for_dim[shuffle([1:numSamples])]
  end
  result
end

end