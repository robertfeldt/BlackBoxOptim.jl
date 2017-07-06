using BlackBoxOptim
using BaseTestNext
using Compat

const NumTestRepetitions = 100

random_candidate(n, lo, hi) =
    BlackBoxOptim.Candidate{Float64}([clamp(randn(), lo, hi) for _ in 1:n])
