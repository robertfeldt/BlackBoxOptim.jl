facts("Population") do

  context("FitPopulation") do
    fs = MinimizingFitnessScheme
    p1 = FitPopulation(fs, 10, 2)
    @fact isa(p1, FitPopulation) => true

    @fact popsize(p1) => 10
    @fact numdims(p1) => 2

    @fact isnafitness(fitness(p1, 1), fs) => true
    @fact isnafitness(fitness(p1, 4), fs) => true

    context("candidates pool") do
      @fact BlackBoxOptim.candi_pool_size(p1) => 0
      candi1 = BlackBoxOptim.acquire_candi(p1, 1)
      @fact BlackBoxOptim.candi_pool_size(p1) => 0
      @fact candi1.index => 1
      @fact isnafitness(candi1.fitness, fs) => true

      candi2 = BlackBoxOptim.acquire_candi(p1, 2)
      @fact BlackBoxOptim.candi_pool_size(p1) => 0
      @fact candi2.index => 2
      @fact isnafitness(candi2.fitness, fs) => true

      BlackBoxOptim.release_candi(p1, candi2)
      @fact BlackBoxOptim.candi_pool_size(p1) => 1

      candi1.fitness = 5.0
      BlackBoxOptim.accept_candi!(p1, candi1)
      @fact BlackBoxOptim.candi_pool_size(p1) => 2
      @fact fitness(p1, 1) => 5.0
    end

    context("append!()") do
      p2 = FitPopulation(fs, 5, 2)
      @fact isa(p2, FitPopulation) => true

      candi22 = BlackBoxOptim.acquire_candi(p2, 2)
      candi22.fitness = 4.0
      BlackBoxOptim.accept_candi!(p2, candi22)

      append!(p1, p2)
      @fact numdims(p1) --> 2
      @fact popsize(p1) --> 15
      @fact p1[11] --> p2[1]
      @fact fitness(p1, 12) --> 4.0

      p3 = FitPopulation(fs, 5, 1)
      @fact_throws DimensionMismatch append!(p1, p3)
    end
  end

end
