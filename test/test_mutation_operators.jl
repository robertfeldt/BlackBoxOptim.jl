facts("Mutation operators") do
  ss = RangePerDimSearchSpace([(-1.0, 1.0), (0.0, 100.0), (-5.0, 0.0)])

  context("SimpleGibbsMutation") do
    gibbs = SimpleGibbsMutation(ss)

    @fact search_space(gibbs) --> ss
    # incorrect dimensions
    @fact_throws BoundsError BlackBoxOptim.apply(gibbs, 2.0, 0, 1)
    @fact_throws BoundsError BlackBoxOptim.apply(gibbs, 2.0, 4, 1)

    for dim in 1:numdims(ss)
      dim_range = range_for_dim(ss, dim)
      for i in 1:NumTestRepetitions
        t = BlackBoxOptim.apply(gibbs, 0.0, dim, 1)
        @fact t --> greater_than_or_equal(dim_range[1])
        @fact t --> less_than_or_equal(dim_range[2])
      end
    end
  end

  context("MutationClock") do
    mc = MutationClock(SimpleGibbsMutation(ss), 0.05)

    mutations_per_param = zeros(numdims(ss))
    NumReps = 2_000
    for i in 1:NumReps
      ref_ind = rand_individual(ss)
      ind = copy(ref_ind)
      BlackBoxOptim.apply!(mc, ind, 1)
      @fact in(ind, ss) --> true
      mutations_per_param[ind .!= ref_ind] += 1
    end
    mut_frequencies = mutations_per_param ./ NumReps
    for p in 1:numdims(ss)
      @fact mut_frequencies[p] --> greater_than(0.03)
      @fact mut_frequencies[p] --> less_than(0.07) # roughly matches the probability
    end
  end

  context("FixedGeneticOperatorsMixture") do
    mx = FixedGeneticOperatorsMixture(GeneticOperator[NoMutation(), MutationClock(SimpleGibbsMutation(ss), 0.05)], [0.3, 0.7])

    ref_ind = rand_individual(ss)
    n_params_mutated = 0
    for i in 1:10000
      ind = copy(ref_ind)
      BlackBoxOptim.apply!(mx, ind, 1)
      @fact in(ind, ss) --> true
      n_params_mutated += ind != ref_ind
    end
    # the number of parameters changed should roughly match the weight of MutationClock multiplied by its mutation probability
    @fact n_params_mutated/numdims(ss) --> greater_than(200)
    @fact n_params_mutated/numdims(ss) --> less_than(500)
  end

  context("FAGeneticOperatorsMixture") do
    mx = FAGeneticOperatorsMixture(GeneticOperator[NoMutation(), MutationClock(SimpleGibbsMutation(ss), 0.05)], pmin = 0.05, pmax = 20.0)

    ref_ind = rand_individual(ss)
    n_params_mutated = 0
    for i in 1:10000
      ind = copy(ref_ind)
      sel_op, tag = BlackBoxOptim.next(mx)
      BlackBoxOptim.apply!(sel_op, ind, 1)
      @fact in(ind, ss) --> true
      is_mutated = ind != ref_ind
      n_params_mutated += is_mutated
      is_improved = rand() < 0.25 && is_mutated
      # fake improvement in 0.25 times mutation clock is applied
      BlackBoxOptim.adjust!(mx, tag, 1, is_improved ? 1.0 : 0.0, 0.0, is_mutated)
    end
    # FA should adjust frequencies to [almost] always apply mutation clock
    # the number of parameters changed should roughly match the MutationClock mutation probability
    @fact n_params_mutated/numdims(ss) --> greater_than(300)
    @fact n_params_mutated/numdims(ss) --> less_than(700)
    @fact frequencies(mx)[1] --> less_than(0.1)
    @fact frequencies(mx)[2] --> greater_than(0.9)
  end

end
