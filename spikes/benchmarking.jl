# We implement Data Profiles for Benchmarking Derivative-Free (i.e. blackbox)
# optimizers as described in the paper:
#
# Jorge J. Moré and Stefan M. Wild (2008), "Benchmarking Derivative-Free
# Optimization Algorithms", Preprint ANL/MCS-P1471-1207, April 2008.
#

# Benchmark a set of optimizers on a set of problems.
# Returns: The data profiles for each of the optimizers on the problems. Each
# data profile uses three performance measures: time, function evaluations and
# iterations (since we always use an optimizer in iterations with one call to
# ask(opt) and one call to tell(opt) for each iteration).
# Both of the first two are relevant and all three can be easily collected
# during optimization. The number of iterations is not so critical to evaluate
# the overall performance of an optimizer, it is more of an internal measure.
# Still it gives an indication of how well an optimizer can be "stepped" and thus
# how flexibly we can use it in more elaborate/complex optimization schemes.
function benchmark(optimizers, problems)

  # We will collect a dict of dicts that map a problem to an opt to the performance
  # trace for that opt on that problem.
  profiles = Dict(Any, Any)

  # For each problem, run all the optimizers while collecting a quality profile.
  # Shuffle the order of the problems just in case they have some dependence.
  for p in shuffle(problems)

    # Create the dict mapping opts to their performance trace on this problem.
    profiles[p] = Dict(Any, Any)

    # Shuffle the order of the optimizers for this problem (in case there is
    # some state/dependence between calls to a problem).
    for opt in shuffle(optimizers)


    end

  end

end

function performance_ratio()
