facts("Latin hypercube sampling") do
  samples = BlackBoxOptim.Utils.latin_hypercube_sampling([0.0], [1.0], 1)
  @fact size(samples, 1) => 1
  @fact 0.0 <= samples[1,1] <= 1.0 => true

  samples = BlackBoxOptim.Utils.latin_hypercube_sampling([0.0], [1.0], 2)
  @fact size(samples, 2) => 2
  sorted = sort(samples, 2)
  @fact 0.0 <= sorted[1,1] <= 0.5 => true
  @fact 0.5 <= sorted[1,2] <= 1.0 => true

  samples = BlackBoxOptim.Utils.latin_hypercube_sampling([0.0, 2.0], [1.0, 3.0], 4)
  @fact size(samples, 1) => 2
  @fact size(samples, 2) => 4
  sorted = sort(samples[1,:],2)
  @fact 0.0 <= sorted[1,1] <= 0.25 => true
  @fact 0.25 <= sorted[1,2] <= 0.5 => true
  @fact 0.5 <= sorted[1,3] <= 0.75 => true
  @fact 0.75 <= sorted[1,4] <= 1.0 => true
  s2 = sort(samples[2,:],2)
  @fact 2.0 <= s2[1,1] <= 2.25 => true
  @fact 2.25 <= s2[1,2] <= 2.5 => true
  @fact 2.5 <= s2[1,3] <= 2.75 => true
  @fact 2.75 <= s2[1,4] <= 3.0 => true
end
