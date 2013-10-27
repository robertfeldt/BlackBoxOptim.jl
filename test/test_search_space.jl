facts("Search space") do
  context("isinspace") do
    for(i in 1:NumTestRepetitions)
      reps = rand(1:100)
      ss1 = symmetric_search_space(reps, (0.0, 1.0))
      ind = rand_individual(ss1)
      for(j in 1:reps)
        @fact (mins(ss1)[j] <= ind[j] <= maxs(ss1)[j])  => true
      end
    end
  end

  context("Symmetric search space with default range") do
    ss1 = symmetric_search_space(1)
    @fact ndims(ss1) => 1
    @fact ranges(ss1) => [(0.0, 1.0)]
    @fact range_for_dim(ss1,1) => (0.0, 1.0)

    for(i in 1:NumTestRepetitions)
      ind = rand_individual(ss1)
      @fact size(ind) => (1, 1)
      @fact isinspace(ind, ss1) => true
    end

    ss3 = symmetric_search_space(3)
    @fact ndims(ss3) => 3
    @fact ranges(ss3) => [(0.0, 1.0), (0.0, 1.0), (0.0, 1.0)]
    @fact range_for_dim(ss3,1) => (0.0, 1.0)
    @fact range_for_dim(ss3,2) => (0.0, 1.0)
    @fact range_for_dim(ss3,3) => (0.0, 1.0)

    for(i in 1:NumTestRepetitions)
      ind = rand_individual(ss3)
      @fact size(ind) => (1, 3)
      @fact isinspace(ind, ss3) => true
    end
  end

  context("SymmetricSearchSpace with given range") do  
    ss1 = symmetric_search_space(1, (-1.0, 1.0))
    @fact ndims(ss1) => 1
    @fact ranges(ss1) => [(-1.0, 1.0)]
    @fact range_for_dim(ss1,1) => (-1.0, 1.0)

    for(i in 1:NumTestRepetitions)
      reps = rand(1:100)
      range = (rand(), rand())
      ss = symmetric_search_space(reps, range)
      @fact ndims(ss) => reps
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
      @fact length(ind) => ndims(ss)
      if !isinspace(ind, ss)
        show(ind)
        show(ss)
      end
      @fact isinspace(ind, ss) => true
    end
  end
end