using BlackBoxOptim: PolynomialMutation, GeneticOperator, MutationOperator, apply

describe("Polynomial mutation genetic operator") do
  @repeat test("mutate float values to new float values within bounds") do
  	low = rand(-10.0:10.0)
  	high = low + rand(0.10:10.0)
  	eta = rand(10.0:250.0)
  	pm = PolynomialMutation(low, high, eta)
  	@check typeof(pm) <: GeneticOperator
  	@check typeof(pm) <: MutationOperator
  	@check typeof(pm) <: PolynomialMutation

  	value = low
  	valueprim = apply(pm, value)
  	@check typeof(valueprim) <: Float64
  	@check low <= valueprim <= high

  	value = high
  	valueprim = apply(pm, value)
  	@check typeof(valueprim) <: Float64
  	@check low <= valueprim <= high

  	for i in 1:100
  		value = low + (high - low) * rand()
  		valueprim = apply(pm, value)
  		@check typeof(valueprim) <: Float64
  		@check low <= valueprim <= high
  	end
  end
end
