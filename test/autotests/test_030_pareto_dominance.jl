using BlackBoxOptim: pareto_dominates, pareto_dominates_hat

describe("Pareto dominance") do
  test("simple 1-objective cases where there is dominance") do
    @check pareto_dominates([0.0], [1.0])
    @check pareto_dominates([-100.0], [100.0])
  end

  test("simple 1-objective cases where there is NO dominance") do
    @check !pareto_dominates([-2.5], [-2.5])
    @check !pareto_dominates([12.32], [12.32])
  end

  test("simple 2-objective cases where there is dominance") do
    @check pareto_dominates([0, 1], [1, 2])
    @check pareto_dominates([-10, -5], [0, 0])
    @check pareto_dominates([0, 1], [0, 2])
    @check pareto_dominates([-1, 2], [0, 2])
  end

  test("simple 2-objective cases where there is NO dominance") do
    @check !pareto_dominates([0, 1], [0, 1])
    @check !pareto_dominates([-1, 1], [-1, 1])

    @check !pareto_dominates([0, 1], [1, 0])
  end

  rand_vector(size = rand(1:107)) = rand(1:103) * randn(size)
  rand_indices(n, numindices) = shuffle(collect(1:n))[1:numindices]

  # Generate a delta-vector for a vector, i.e. a vector of changes in a number of nonzero positions.
  delta_vector_for(u, numNonzero) = begin
    size = length(u)
    delta = abs(rand_vector(size))
    delta[rand_indices(size, size - min(numNonzero, size))] = 0.0
    delta
  end

  @repeat test("when the vectors are equal, i.e. NO dominance") do
    u = rand_vector()
    @check !pareto_dominates(u, u)
    @check pareto_dominates_hat(u, u) == 0
  end

  @repeat test("when the vectors have one larger and one smaller value, i.e. NO dominance") do
    u = rand_vector(rand(2:107))
    v = copy(u)
    idxs = rand_indices(length(u), 2)
    v[idxs[1]] += rand()
    v[idxs[2]] -= rand()
    @check !pareto_dominates(u, v)
    @check pareto_dominates_hat(u, v) == 0
  end

  @repeat test("when there is one dominating") do
    u = rand_vector()
    v = copy(u) .+ abs(delta_vector_for(u, rand(1:length(u))))
    @check pareto_dominates(u, v)
    @check !pareto_dominates(v, u)
    @check pareto_dominates_hat(u, v) == -1
    @check pareto_dominates_hat(v, u) ==  1
  end

  @repeat test("for single float fitness values") do
    fval1 = randn()
    f1 = makefitness(fval1)
    fval2 = fval1 + rand(0.0001:113.0)
    f2 = makefitness(fval2)
    @check pareto_dominates(f1, f2)
    @check !pareto_dominates(f2, f1)
    @check !pareto_dominates(f1, f1)
    @check !pareto_dominates(f2, f2)

    fval3 = fval1 - rand(0.00001:113.0)
    f3 = makefitness(fval3)
    @check pareto_dominates(f3, f1)
    @check pareto_dominates(f3, f2)
    @check !pareto_dominates(f1, f3)
    @check !pareto_dominates(f2, f3)
    @check !pareto_dominates(f3, f3)
  end
end