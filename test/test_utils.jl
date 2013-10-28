facts("Utilities") do
  context("Latin hypercube sampling") do
    samples = GlobalOptim.Utils.latin_hypercube_sampling([0.0], [1.0], 1)
    @fact size(samples, 1) => 1
    @fact 0.0 <= samples[1,1] <= 1.0 => true

    samples = GlobalOptim.Utils.latin_hypercube_sampling([0.0], [1.0], 2)
    @fact size(samples, 1) => 2
    sorted = sort(samples[:,1])
    @fact 0.0 <= sorted[1,1] <= 0.5 => true
    @fact 0.5 <= sorted[2,1] <= 1.0 => true

    samples = GlobalOptim.Utils.latin_hypercube_sampling([0.0 2.0], [1.0 3.0], 2)
    @fact size(samples, 1) => 2
    @fact size(samples, 2) => 2
    sorted = sort(samples, 1)
    @fact 0.0 <= sorted[1,1] <= 0.5 => true
    @fact 0.5 <= sorted[2,1] <= 1.0 => true
    s2 = sort(samples[:,2])
    @fact 2.0 <= sorted[1,2] <= 2.5 => true
    @fact 2.5 <= sorted[2,2] <= 3.0 => true
  end
end