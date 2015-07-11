facts("Population") do

  context("FitPopulation") do
    fs = ScalarFitness{true}()
    p1 = FitPopulation(fs, 10, 2)
    @fact isa(p1, FitPopulation) => true

    @fact popsize(p1) => 10
    @fact numdims(p1) => 2

    @fact isnafitness(fitness(p1, 1), fs) => true
    @fact isnafitness(fitness(p1, 4), fs) => true

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

end
