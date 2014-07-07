# Direct Search with Probabilistic Descent as described in Gratton2014:
#
#  S. Gratton, C. W. Royer, L. N. Vicente, and Z. Zhang, Direct search 
#  based on probabilistic descent, preprint 14-11, Dept. Mathematics, Univ. Coimbra.
#  http://www.mat.uc.pt/~lnv/papers/ds-random.pdf
#

# Generate num random vectors on the n-dimensional, unit (hyper)sphere.
# This is the Muller-Marsaglia method as described on the page:
#   http://mathworld.wolfram.com/HyperspherePointPicking.html
function sample_unit_hypersphere(n, num = 1)
  X = randn(n, num)
  sqrootsums = 1 ./ sqrt(sum( X.^2, 1 ))
  broadcast(*, sqrootsums, X)
end

type RandomDirectionGen <: DirectionGenerator
  numDimensions::Int64
  numDirections::Int64
end

function directions_for_k(rdg::RandomDirectionGen, k)
  sample_unit_hypersphere(rdg.numDimensions, rdg.numDirections)
end

# A MirroredRandomDirectionGen generates half of the directions randomly and then
# mirrors by negating them.
type MirroredRandomDirectionGen <: DirectionGenerator
  numDimensions::Int64
  numDirections::Int64
end

function directions_for_k(rdg::RandomDirectionGen, k)
  sample_unit_hypersphere(rdg.numDimensions, rdg.numDirections)
end
