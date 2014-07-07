facts("Random sampling on unit, n-dimensional sphere") do

  uv1 = BlackBoxOptim.sample_unit_hypersphere(1)
  @fact size(uv1, 1) => 1

  uv2 = BlackBoxOptim.sample_unit_hypersphere(2)
  @fact size(uv2, 1) => 2

  uv3 = BlackBoxOptim.sample_unit_hypersphere(3, 4)
  @fact size(uv3, 1) => 3
  @fact size(uv3, 2) => 4

  for(i in 1:100)
    n = rand(1:100)
    num = rand(2:10)
    u = BlackBoxOptim.sample_unit_hypersphere(n, num)
    @fact size(u, 1) => n
    @fact size(u, 2) => num
    for(j in 1:num)
      @fact isapprox( norm(u[:,j]), 1.0 ) => true
    end
  end

end