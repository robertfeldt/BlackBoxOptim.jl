"""
Simple random `IndividualsSelector`.

The probabilties of all candidates are equal.
"""
struct SimpleSelector <: IndividualsSelector
end

select(::SimpleSelector, population, numSamples::Int) =
    rand_indexes(1:popsize(population), numSamples)
