facts("Search space") do
  context("isinspace") do
    for(i in 1:NumTestRepetitions)
      reps = rand(1:10)
      ss1 = symmetric_search_space(reps, (0.0, 1.0))
      ind = rand_individual(ss1)
      for(j in 1:reps)
        @fact (mins(ss1)[j] <= ind[j] <= maxs(ss1)[j])  => true
      end
    end
  end

  context("Symmetric search space with default range") do
    ss1 = symmetric_search_space(1)
    @fact numdims(ss1) => 1
    @fact ranges(ss1) => [(0.0, 1.0)]
    @fact range_for_dim(ss1,1) => (0.0, 1.0)

    for(i in 1:NumTestRepetitions)
      ind = rand_individual(ss1)
      @fact size(ind) => (1,)
      @fact isinspace(ind, ss1) => true
    end

    ss3 = symmetric_search_space(3)
    @fact numdims(ss3) => 3
    @fact ranges(ss3) => [(0.0, 1.0), (0.0, 1.0), (0.0, 1.0)]
    @fact range_for_dim(ss3,1) => (0.0, 1.0)
    @fact range_for_dim(ss3,2) => (0.0, 1.0)
    @fact range_for_dim(ss3,3) => (0.0, 1.0)

    for(i in 1:NumTestRepetitions)
      ind = rand_individual(ss3)
      @fact size(ind) => (3,)
      @fact isinspace(ind, ss3) => true
    end
  end

  context("SymmetricSearchSpace with given range") do
    ss1 = symmetric_search_space(1, (-1.0, 1.0))
    @fact numdims(ss1) => 1
    @fact ranges(ss1) => [(-1.0, 1.0)]
    @fact range_for_dim(ss1,1) => (-1.0, 1.0)

    for(i in 1:NumTestRepetitions)
      reps = rand(1:100)
      range = (rand(), rand())
      ss = symmetric_search_space(reps, range)
      @fact numdims(ss) => reps
      @fact all([(dr == range) for dr in ranges(ss)]) => true
    end
  end

  context("rand_individual is within the search space") do
    for(i in 1:NumTestRepetitions)
      reps = rand(1:100)
      mm = sort(rand(2,1), 1)
      range = (mm[1], mm[2])
      ss = symmetric_search_space(reps, range)
      ind = rand_individual(ss)
      @fact length(ind) => numdims(ss)
      @fact isinspace(ind, ss) => true
    end
  end

  context("rand_individuals creates many individuals and all are within the search space") do
    for(i in 1:NumTestRepetitions)
      reps = rand(1:10)
      mm = sort(rand(2,1), 1)
      range = (mm[1], mm[2])
      ss = symmetric_search_space(reps, range)
      numinds = rand(1:10)
      inds = rand_individuals(ss, numinds)
      @fact size(inds,1) => numdims(ss)
      @fact size(inds,2) => numinds
      for(j in 1:numinds)
        @fact isinspace(inds[:,j], ss) => true
      end
    end
  end

  context("rand_individuals correctly handles column-wise generation in assymetric search spaces") do
    for(i in 1:NumTestRepetitions÷10)
      numdims = rand(1:13)
      minbounds = rand(numdims)
      ds = rand(1:10, numdims) .* rand(numdims)
      maxbounds = minbounds .+ ds
      parambounds = collect(zip(minbounds, maxbounds))
      ss = RangePerDimSearchSpace(parambounds)
      @fact mins(ss) => minbounds
      @fact maxs(ss) => maxbounds
      @fact round(deltas(ss), 6) => round(ds, 6)

      # Now generate 100 individuals and make sure they are all within bounds
      inds = rand_individuals(ss, 100)
      @fact size(inds, 2) => 100
      for i in 1:size(inds, 2)
        for d in 1:numdims
          @fact (minbounds[d] <= inds[d,i] <= maxbounds[d]) => true
        end
      end
    end
  end

  context("RangePerDimSearchSpace") do
    ss = RangePerDimSearchSpace([(0.0, 1.0)])
    @fact mins(ss) => [0.0]
    @fact maxs(ss) => [1.0]
    @fact deltas(ss) => [1.0]

    ss = RangePerDimSearchSpace([(0.0, 1.0), (0.5, 10.0)])
    @fact mins(ss) => [0.0, 0.5]
    @fact maxs(ss) => [1.0, 10.0]
    @fact deltas(ss) => [1.0, 9.5]
  end

  context("rand_individuals_lhs samples in LHS intervals") do
    ss = RangePerDimSearchSpace([(0.0, 1.0), (2.0, 3.0), (4.0, 5.0)])

    inds = rand_individuals_lhs(ss, 2)
    @fact size(inds, 1) => 3
    @fact size(inds, 2) => 2

    sorted = sort(inds, 2) # Sort per row => in their ordered intervals
    @fact (0.0 <= sorted[1,1] <= 0.5) => true
    @fact (0.5 <= sorted[1,2] <= 1.0) => true

    @fact (2.0 <= sorted[2,1] <= 2.5) => true
    @fact (2.5 <= sorted[2,2] <= 3.0) => true

    @fact (4.0 <= sorted[3,1] <= 4.5) => true
    @fact (4.5 <= sorted[3,2] <= 5.0) => true
  end

  context("feasible finds feasible points in the search space") do
    ss = RangePerDimSearchSpace([(0.0, 1.0), (2.0, 3.0), (4.0, 5.0)])

    # We use the double transpose below to ensure the actual and expected
    # values have the same type (matrices, not vectors).
    @fact BlackBoxOptim.feasible([1.1, 2.0, 4.0], ss) => [1.0, 2.0, 4.0]
    @fact BlackBoxOptim.feasible([1.1, 3.0, 4.0], ss) => [1.0, 3.0, 4.0]
    @fact BlackBoxOptim.feasible([1.1, 2.0, 5.0], ss) => [1.0, 2.0, 5.0]
    @fact BlackBoxOptim.feasible([1.1, 3.0, 5.0], ss) => [1.0, 3.0, 5.0]

    @fact BlackBoxOptim.feasible([-0.1, 2.0, 4.0], ss) => [0.0, 2.0, 4.0]
    @fact BlackBoxOptim.feasible([-0.1, 3.0, 4.0], ss) => [0.0, 3.0, 4.0]
    @fact BlackBoxOptim.feasible([-0.1, 2.0, 5.0], ss) => [0.0, 2.0, 5.0]
    @fact BlackBoxOptim.feasible([-0.1, 3.0, 5.0], ss) => [0.0, 3.0, 5.0]

    @fact BlackBoxOptim.feasible([0.0, 1.9, 4.0], ss) => [0.0, 2.0, 4.0]
    @fact BlackBoxOptim.feasible([0.0, 1.9, 4.0], ss) => [0.0, 2.0, 4.0]
    @fact BlackBoxOptim.feasible([1.0, 1.9, 5.0], ss) => [1.0, 2.0, 5.0]
    @fact BlackBoxOptim.feasible([1.0, 1.9, 5.0], ss) => [1.0, 2.0, 5.0]

    @fact BlackBoxOptim.feasible([0.0, 3.3, 4.0], ss) => [0.0, 3.0, 4.0]
    @fact BlackBoxOptim.feasible([0.0, 3.2, 4.0], ss) => [0.0, 3.0, 4.0]
    @fact BlackBoxOptim.feasible([1.0, 3.1, 5.0], ss) => [1.0, 3.0, 5.0]
    @fact BlackBoxOptim.feasible([1.0, 3.9, 5.0], ss) => [1.0, 3.0, 5.0]

    @fact BlackBoxOptim.feasible([-0.4, 3.3, 14.5], ss) => [0.0, 3.0, 5.0]
  end

  context("diameters") do
    ss = RangePerDimSearchSpace([(0.0, 1.0), (2.0, 3.0), (4.0, 5.0)])
    diams = diameters(ss)

    @fact length(diams) => 3
    @fact diams[:] => [1.0, 1.0, 1.0]
  end

end
