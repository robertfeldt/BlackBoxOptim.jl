facts("Random search") do

context("ask") do
  for(i in 1:NumTestRepetitions)

    dims = rand(1:100)
    range_per_dim = (0.0, 1.0)
    ss = [range_per_dim for i in 1:dims]

    opt = random_searcher(ss)

    res = GlobalOptim.ask(opt)

    @fact length(res) => 1
    trial, trialIndex = res[1]

    @fact ndims(trial) => dims
    @fact trialIndex => Nothing
    @fact is_within_bounds(trial, ss) => true

  end
end

end