facts("Random search") do

context("ask") do
  for(i in 1:NumTestRepetitions)
    dims = rand(1:20)
    min = rand(1:123)
    range = (min * rand(), min + rand() * min)

    ss = BlackBoxOptim.symmetric_search_space(dims, range)
    opt = BlackBoxOptim.random_search(ss)

    res1 = BlackBoxOptim.ask(opt)

    @fact length(res1) => 1
    candidate, index = res1[1]
    @fact size(candidate,2) => dims
    @fact isinspace(candidate, ss) => true

    # Fake fitness
    better = BlackBoxOptim.tell(opt, [(candidate, index, 1.0)])

    @fact better => 1
    @fact opt.best => candidate
    @fact opt.best_fitness => 1.0

    # Get one more and fake that it has better fitness.
    res2 = BlackBoxOptim.ask(opt)
    @fact length(res2) => 1
    candidate2, index2 = res2[1]
    @fact size(candidate2,2) => dims
    @fact isinspace(candidate2, ss) => true
    better2 = BlackBoxOptim.tell(opt, [(candidate, index, 0.5)])
  end
end

end