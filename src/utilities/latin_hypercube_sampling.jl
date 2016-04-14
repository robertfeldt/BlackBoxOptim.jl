"""
    latin_hypercube_sampling(mins, maxs, numSamples)

    Randomly sample `numSamples` values from the parallelogram defined
    by `mins` and `maxs` using the Latin hypercube algorithm.
"""
function latin_hypercube_sampling{T}(mins::AbstractVector{T}, maxs::AbstractVector{T}, numSamples::Integer)
    dims = length(mins)
    result = zeros(T, numSamples, dims)
    @inbounds for i in 1:dims
        interval_len = (maxs[i] - mins[i]) / numSamples
        result[:,i] = shuffle!(linspace(mins[i], maxs[i] - interval_len, numSamples) +
                               interval_len*rand(numSamples))
    end
    return result'
end
