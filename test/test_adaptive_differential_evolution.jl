NumTestRepetitions = 100

facts("Adaptive differential evolution optimizer") do

ade = adaptive_de_rand_1_bin()

context("sample_f") do
  for(i in 1:NumTestRepetitions)
    @fact 0.0 <= GlobalOptim.sample_f(ade) <= 1.0 => true
  end
end

context("ask") do
  for(i in 1:NumTestRepetitions)
    res = GlobalOptim.ask(ade)

    @fact length(res) => 2
    trial, trialIndex = res[1]
    target, targetIndex = res[2]

    @fact ndims(trial) => 2
    @fact 1 <= trialIndex <= length(ade.population) => true
    @fact isinspace(trial, ade.search_space) => true

    @fact ndims(target) => 2
    @fact 1 <= targetIndex <= length(ade.population) => true
    @fact isinspace(target, ade.search_space) => true

    @fact trialIndex == targetIndex => true
  end
end

end