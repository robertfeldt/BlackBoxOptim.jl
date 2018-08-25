"""
    latin_hypercube_sampling(mins, maxs, n)

Randomly sample `n` vectors from the parallelogram defined
by `mins` and `maxs` using the Latin hypercube algorithm.

Returns `dims`×`n` matrix.
"""
function latin_hypercube_sampling(mins::AbstractVector{T},
                                  maxs::AbstractVector{T},
                                  n::Integer) where T<:Number
    length(mins) == length(maxs) ||
        throw(DimensionMismatch("mins and maxs should have the same length"))
    all(xy -> xy[1] <= xy[2], zip(mins, maxs)) ||
        throw(ArgumentError("mins[i] should not exceed maxs[i]"))
    dims = length(mins)
    result = zeros(T, n, dims)
    @inbounds for i in 1:dims
        interval_len = (maxs[i] - mins[i]) / n
        result[:,i] = shuffle!(range(mins[i], stop=maxs[i] - interval_len, length=n) +
                               interval_len*rand(n))
    end
    return result'
end
