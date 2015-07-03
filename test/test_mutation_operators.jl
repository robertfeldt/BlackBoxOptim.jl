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

    n_params_mutated = 0
    for i in 1:10000
      ref_ind = rand_individual(ss)
      ind = copy(ref_ind)
      BlackBoxOptim.apply!(mc, ind)
      @fact isinspace(ind, ss) => true
      n_params_mutated += sum(ind .!= ref_ind)
    end
    @fact 300 <= n_params_mutated/numdims(ss) <= 700 => true # roughly matches the probability
  end

  context("MutationMixture") do
    mx = MutationMixture([NoMutation(), MutationClock(SimpleGibbsMutation(ss), 0.05)], [0.3, 0.7])

    ref_ind = rand_individual(ss)
    n_params_mutated = 0
    for i in 1:10000
      ind = copy(ref_ind)
      BlackBoxOptim.apply!(mx, ind)
      @fact isinspace(ind, ss) => true
      n_params_mutated += (any(ind .!= ref_ind))
    end
    # the number of parameters changed should roughly match the weight of MutationClock multiplied by its mutation probability
    @fact 200 < n_params_mutated/numdims(ss) < 500 => true
  end

end
