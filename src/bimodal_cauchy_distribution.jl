using Distributions

# FIXME implement actual distribution using Distributions.MixtureModel
"""
A mixture of 2 Cauchy distributions.
Random values are further constrained to `[0.0, 1.0]` range either by
truncating the initial unconstrained value or by generating new random value
until it fits the range.
"""
struct BimodalCauchy
    a::Cauchy
    b::Cauchy
    mix_prob::Float64
    # When sampling it is common to truncate in either or both ends.
    clampBelow0::Bool
    clampAbove1::Bool

    BimodalCauchy(location1, scale1, location2, scale2;
                  mix_prob::Number=0.5, clampBelow0::Bool=true, clampAbove1::Bool=true) =
        new(Cauchy(location1, scale1), Cauchy(location2, scale2),
            mix_prob, clampBelow0, clampAbove1)
end

function Random.rand(distr::BimodalCauchy)
    while true
        value = rand() < distr.mix_prob ? rand(distr.a) : rand(distr.b)
        if value >= 1.0
            distr.clampAbove1 && return 1.0
        elseif value <= 0.0
            distr.clampBelow0 && return 0.0
        else
            return value
        end
    end
end

function Random.rand!(distr::BimodalCauchy, A::AbstractArray)
    for i in eachindex(A)
        @inbounds A[i] = rand(distr)
    end
    return A
end
