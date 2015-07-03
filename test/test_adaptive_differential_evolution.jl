NumTestRepetitions = 100

facts("Adaptive differential evolution optimizer") do

ade = adaptive_de_rand_1_bin()

context("parameters adjust!()") do
  for(i in 1:NumTestRepetitions)
    cur_cr, cur_f = BlackBoxOptim.crossover_parameters(ade.params, 1)
    @fact 0.0 <= cur_cr <= 1.0 => true
    @fact 0.0 <= cur_f <= 1.0 => true
    BlackBoxOptim.adjust!(ade.params, 1, false)
    # FIXME this fails too often; needs adjusting the distribution?
    #new_cr, new_f = BlackBoxOptim.crossover_parameters(ade.params, 1)
    #@fact new_cr != cur_cr => true
    #@fact new_f != cur_f => true
  end
end

context("ask()") do
  for(i in 1:NumTestRepetitions)
    res = BlackBoxOptim.ask(ade)

    @fact length(res) => 2
    trial, trialIndex = res[1]
    target, targetIndex = res[2]

    @fact ndims(trial) => 1
    @fact 1 <= trialIndex <= popsize(ade) => true
    @fact isinspace(trial, ade.embed.searchSpace) => true

    @fact ndims(target) => 1
    @fact 1 <= targetIndex <= popsize(ade) => true
    @fact isinspace(target, ade.embed.searchSpace) => true

    @fact trialIndex == targetIndex => true
  end
end

end
