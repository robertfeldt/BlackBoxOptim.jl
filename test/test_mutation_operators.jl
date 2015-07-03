facts("Mutation operators") do
  ss = RangePerDimSearchSpace([(-1.0, 1.0), (0.0, 100.0), (-5.0, 0.0)])

  context("SimpleGibbsMutation") do
    gibbs = SimpleGibbsMutation(ss)

    @fact search_space(gibbs) => ss
    # incorrect dimensions
    @fact_throws BoundsError BlackBoxOptim.apply(gibbs, 2.0, 0)
    @fact_throws BoundsError BlackBoxOptim.apply(gibbs, 2.0, 4)

    for dim in 1:numdims(ss)
      dim_range = range_for_dim(ss, dim)
      for i in 1:NumTestRepetitions
        t = BlackBoxOptim.apply(gibbs, 0.0, dim)
        @fact dim_range[1] <= t <= dim_range[2] => true
      end
    end
  end

  context("MutationClock") do
    mc = MutationClock(SimpleGibbsMutation(ss), 0.05)

    for i in 1:NumTestRepetitions
      ref_ind = rand_individual(ss)
      ind = copy(ref_ind)
      BlackBoxOptim.apply!(mc, ind)
      @fact isinspace(ind, ss) => true
      @fact any(ind .!= ref_ind) => true
      @fact all(ind .!= ref_ind) => false # might be true with very low probability
    end
  end

end
