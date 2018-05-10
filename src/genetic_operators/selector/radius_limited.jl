"""
`IndividualsSelector` that implements a "trivial geography" similar to Spector and Kline (2006)
by first sampling an individual randomly and then selecting additional
individuals for the same tournament within a certain deme of limited size (`radius`)
for the sub-sequent individuals in the population.

The version we implement here is from:
    I. Harvey, "The Microbial Genetic Algorithm", in Advances in Artificial Life
    Darwin Meets von Neumann, Springer, 2011.
The original paper is:
    Spector, L., and J. Klein. 2005. Trivial Geography in Genetic Programming.
    In Genetic Programming Theory and Practice III, edited by T. Yu, R.L. Riolo,
    and B. Worzel, pp. 109-124. Boston, MA: Kluwer Academic Publishers.
    http://faculty.hampshire.edu/lspector/pubs/trivial-geography-toappear.pdf
"""
mutable struct RadiusLimitedSelector <: IndividualsSelector
    radius::Int
end

function select(sel::RadiusLimitedSelector, population, numSamples::Int)
    # The radius must be at least as big as the number of samples + 2 so that
    # there is something to sample from.
    radius = max(sel.radius, numSamples+2)
    psize = popsize(population)
    deme_start = rand(1:psize)
    ixs = rand_indexes(deme_start:(deme_start+radius-1), numSamples)
    # Ensure they are not out of bounds by wrapping over at the end.
    ixs .= mod1.(ixs, psize)
end
