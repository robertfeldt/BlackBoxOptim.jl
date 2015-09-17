facts("Random search") do

context("ask()") do
  for(i in 1:NumTestRepetitions)
    dims = rand(1:20)
    min = rand(1:123)
    range = (min * rand(), min + rand() * min)

    ss = BlackBoxOptim.symmetric_search_space(dims, range)
    opt = BlackBoxOptim.random_search(ss)

    res1 = BlackBoxOptim.ask(opt)

    @fact length(res1) --> 1
    candidate = res1[1]
    @fact length(candidate.params) --> dims
    @fact in(candidate.params, ss) --> true

    # Fake fitness
    candidate.fitness = 1.0
    better = BlackBoxOptim.tell!(opt, [candidate])

    @fact better --> 1
    @fact opt.best --> candidate.params
    @fact opt.best_fitness --> 1.0

    # Get one more and fake that it has better fitness.
    res2 = BlackBoxOptim.ask(opt)
    @fact length(res2) --> 1
    candidate2 = res2[1]
    @fact length(candidate2.params) --> dims
    @fact in(candidate2.params, ss) --> true
    candidate2.fitness = 0.5
    better2 = BlackBoxOptim.tell!(opt, [candidate2])
  end
end

end
