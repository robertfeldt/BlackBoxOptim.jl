# Archives work with candidate solutions. All candidate solutions can be mapped to
# a behavioral descriptor (aka a fitness vector) which is used by the archive to judge 
# what to maintain in the archive.
abstract Candidate

behavior(c::Candidate) = c.behavior

# The most typical candidate is a vector of numbers to represent the genotype
# and a vector of numbers to represent the behavior of the phenotype.
type VectorCandidate{ET <: Real, BT <: Real} <: Candidate
  genotype::Vector{ET}
  behavior::Union(Nothing, Vector{BT})
  VectorCandidate{ET <: Real, BT <: Real}
end

# An epsilon-box archive saves only the solutions that are not epsilon-box
# dominated by any other solutions in the archive. It also counts the number
# of candidate solutions that have been added and how many epsilon-box progresses
# that have been made.
type EpsilonBoxArchive <: Archive
end

# An EpsilonBoxArchive saves candidates in the archive. Each candidate must implement
# a fitness_vector function that returns a fitness vector, i.e. float values that 
# represents the objective values of the candidate.
fitness_vector()