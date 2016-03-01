"""
  Simple random `IndividualsSelector`.

  The probabilties of all candidates are equal.
"""
immutable SimpleSelector <: IndividualsSelector
end

function select(::SimpleSelector, population, numSamples::Int)
  sample(1:popsize(population), numSamples; replace = false)
end
