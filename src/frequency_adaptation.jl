# A FrequencyAdapter adapts the frequencies with which a set of 
# values/strategies/methods should be applied/tried in an optimization problem.
# It is based on the Adaptive Coordinate Frequencies scheme described in:
#
# T. Glasmachers and U. Dogan, "Accelerated Coordinate Descent with 
# Adaptive Coordinate Frequencies", 2013.
#
# but generalized so that it can support more than the adaptation of only
# the coordinates in a Coordinate Descent scheme. The things that are being
# adapted are identified by integers in 1:n, with n being the main parameter.
type FrequencyAdapter
  n::Integer
  eta::Float64
  c::Float64
  pmin::Float64
  pmax::Float64

  p::Vector{Float64}    # Current pi values.
  psum::Float64         # Current psum.
  a::Vector{Float64}    # Current ai values.
  deltahat::Float64     # Running average of the progress values.

  block::Vector{Int}  # Current block of selected methods.

  numupdates::Int     # Number of times we have been updated.
  min_updates::Int    # Number of updates until we start adapting frequencies.

  FrequencyAdapter(n, c = 0.2, pmin = 0.05, pmax = 20) = begin
    eta = 1/n
    p = ones(Float64, n)/n
    new(n, eta, c, pmin, pmax, 
      p, 1.0, 
      zeros(Float64, n),
      0.0,
      shuffle(collect(1:n)),
      0, n)
  end
end

# Return the index of the next method that should be used. The base 
# FrequencyAdapter first creates a block of randomly shuffled methods that should
# be applied and then selected the next one from it. If there is a current
# block which is non-empty, use that, if not create a new block. However,
# the first block is always random shuffles of all the methods since we need
# to learn about their effectiveness.
function next(fa::FrequencyAdapter)
  if length(fa.block) < 1
    create_new_block!(fa)
  end
  #println("Taking last from block = ", fa.block)
  return pop!(fa.block)
end

# Create a new block.
function create_new_block!(fa::FrequencyAdapter)
  block = Int[]
  #print("Creating new block, psum = $(fa.psum), a = ", fa.a, ", p = ", fa.p)
  for(i in 1:fa.n)
    fa.a[i] += (fa.n * fa.p[i] / fa.psum)
    num_ai = convert(Int, floor(fa.a[i]))
    if num_ai > 0
      block = vcat(block, ones(Int, num_ai)*i)
      fa.a[i] -= num_ai
    end
  end
  # Due to rounding errors the block is sometimes empty so we select the one 
  # with largest a value.
  if length(block) < 1
    #println("Empty block created. Rectifying., std(p) = ", std(fa.p))
    #print("psum = $(fa.psum), a = ", fa.a, ", p = ", fa.p)
    i = indmax(fa.a)
    fa.a[i] = 0
    block = [i]
  end
  #println("Created new block = ", block)
  fa.block = shuffle(block)
end

# Update the internal model of progress and success rate of each method based
# on the latest progress value of one method. Progress values should be larger
# the larger the progress/improvement was.
function update!(fa::FrequencyAdapter, methodIndex, progress)
  # If we already have collected a few samples of progress rates we can update
  # the pi. This is the common case.
  if fa.numupdates >= fa.min_updates
    unconstrained_pnew = fa.p[methodIndex] * exp(fa.c * (progress / fa.deltahat - 1))
    pnew = min(fa.pmax, max(fa.pmin, unconstrained_pnew))
    #print("i = ", methodIndex, ", pi = ", fa.p[methodIndex], ", pnew = ", pnew, ", psum = ", fa.psum)
    fa.psum += (pnew - fa.p[methodIndex])
    fa.p[methodIndex] = pnew
    #print(", new psum = ", fa.psum, ", deltahat = ", fa.deltahat)
    fa.deltahat = (1 - fa.eta) * fa.deltahat + fa.eta * progress
    #println(", new deltahat = ", fa.deltahat)
  else
    # Until we have collected at least min_updates we just sum the progress 
    # values so we can later calculate their average.
    fa.deltahat += progress
  end
  fa.numupdates += 1
  if fa.numupdates == fa.min_updates
    #print("Setting deltahat, was = ", fa.deltahat, " (min_updates = $(fa.min_updates))")
    fa.deltahat = fa.deltahat / fa.min_updates
    #println(", now = ", fa.deltahat)
  end
end