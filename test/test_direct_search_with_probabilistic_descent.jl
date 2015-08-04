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

facts("Random direction generator") do

  rdg1 = BlackBoxOptim.RandomDirectionGen(2, 3)
  ds1 = BlackBoxOptim.directions_for_k(rdg1, 1)
  @fact size(ds1) => (2, 3)

  rdg2 = BlackBoxOptim.RandomDirectionGen(10, 17)
  ds2 = BlackBoxOptim.directions_for_k(rdg2, 1)
  @fact size(ds2) => (10, 17)

end

facts("Mirrored random direction generator") do

  mrdg1 = BlackBoxOptim.MirroredRandomDirectionGen(2, 4)
  ds1 = BlackBoxOptim.directions_for_k(mrdg1, 1)
  @fact size(ds1) => (2, 4)
  @fact ds1[:,3] => -ds1[:,1]
  @fact ds1[:,4] => -ds1[:,2]

  mrdg2 = BlackBoxOptim.MirroredRandomDirectionGen(10, 6)
  ds2 = BlackBoxOptim.directions_for_k(mrdg2, 1)
  @fact size(ds2) => (10, 6)
  @fact ds2[:,4] => -ds2[:,1]
  @fact ds2[:,5] => -ds2[:,2]
  @fact ds2[:,6] => -ds2[:,3]

  # Must be even number of directions
  @fact_throws BlackBoxOptim.MirroredRandomDirectionGen(10, 1)
  @fact_throws BlackBoxOptim.MirroredRandomDirectionGen(10, 3)
  @fact_throws BlackBoxOptim.MirroredRandomDirectionGen(10, 7)

end
