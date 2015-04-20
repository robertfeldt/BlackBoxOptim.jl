using BlackBoxOptim: MutationClock, apply, PolynomialMutation

count_mutations(origvector, mutatedvector) = sum(origvector .!= mutatedvector)

describe("count_mutations") do
  test("counts changes") do
    @check count_mutations([0,0], [0,0]) == 0
    @check count_mutations([0,0], [1,0]) == 1
    @check count_mutations([0,0], [0,1]) == 1
    @check count_mutations([0,0], [1,1]) == 2
    @check count_mutations([0,1], [1,1]) == 1
    @check count_mutations([1,1], [1,1]) == 0

    @check count_mutations([0.0, 0.0], [1.8, 0.0]) == 1
    @check count_mutations([0.0, 11.32], [1.8, 0.0]) == 2
  end
end

describe("Mutation Clock genetic operator") do
  @repeat test("expected number of mutations are done (if passed vectors of same size)") do
    probmutate = rand(0.000001:0.001:1.0)
    pmut = PolynomialMutation(0.0, 1.0)
    mc = MutationClock(pmut, probmutate)
    size = rand(1:100)

    # Now apply mutation clock operator to a large number of vectors and count the total
    # number of mutations
    nmutations = 0
    nreps = 100
    for i in 1:nreps
      vector = rand(size)
      nmutations += count_mutations(vector, apply(mc, copy(vector)))
    end

    # Not sure why the interval must be so wide. Seems we always get lower than the
    # expected number of mutations.
    nexpected = probmutate * size * nreps
    @check ifloor(0.50 * nexpected) <= nmutations <= iceil(1.50 * nexpected)
  end
end