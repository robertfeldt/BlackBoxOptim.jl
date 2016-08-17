"""
    A point in the problem's search space with the
    corresponding fitness value.

  `F` is the original problem's fitness type
"""
abstract FitIndividual{F}

fitness_type{F}(indi::FitIndividual{F}) = F

"""
    Get the problem parameters (a point in the search space) of the individual.
"""
params(indi::FitIndividual) = indi.params

"""
    fitness(indi::FitIndividual)

    Gets the fitness of the individual.
"""
fitness(indi::FitIndividual) = indi.fitness

"""
  Candidate solution to the problem.

  Candidate can be either a member of the population (`index` > 0) or
  a standalone solution (`index` == -1).
  Can carry additional information, like the `tag` or the genetic operator applied (`extra`).
"""
type Candidate{F} <: FitIndividual{F}
    params::Individual
    index::Int          # index of individual in the population, -1 if unassigned
    fitness::F          # fitness

    extra::Any          # extra information
    tag::Int            # additional information set by the genetic operator

    Candidate(params::Individual, index::Int = -1,
             fitness::F = NaN,
             extra::Any = Void,
             tag::Int = 0) =
        new(params, index, fitness, extra, tag)

    @compat (::Type{Candidate}){F}(
            params::Individual, index::Int = -1,
            fitness::F = NaN,
            extra::Any = Void,
            tag::Int = 0) =
        new{F}(params, index, fitness, extra, tag)
end

index(cand::Candidate) = cand.index
tag(cand::Candidate) = cand.tag

Base.copy(c::Candidate) = Candidate(copy(c.params), c.index, c.fitness, c.extra, c.tag)

function Base.copy!{F}(c::Candidate{F}, o::Candidate{F})
    copy!(c.params, o.params)
    c.index = o.index
    c.fitness = o.fitness # FIXME if vector?
    c.extra = o.extra
    c.tag = o.tag
    return c
end
