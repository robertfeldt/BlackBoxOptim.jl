facts("Bimodal Cauchy Distributions") do

  context("sample bimodal cauchy with truncation on one side") do
    bc = BlackBoxOptim.bimodal_cauchy(0.65, 0.1, 1.0, 0.1)

    for i in 1:100
      v = BlackBoxOptim.sample_bimodal_cauchy(bc; truncateAbove1 = true, truncateBelow0 = false)

      @fact v --> less_than_or_equal(1.0)
      @fact v --> greater_than(0.0) # Very unlikely to be 0.0!!?
    end
  end

  context("sample bimodal cauchy with truncation on both sides") do
    bc = BlackBoxOptim.bimodal_cauchy(0.1, 0.1, 0.95, 0.1)

    for i in 1:100
      v = BlackBoxOptim.sample_bimodal_cauchy(bc; truncateAbove1 = true, truncateBelow0 = true)

      @fact v --> less_than_or_equal(1.0)
      @fact v --> greater_than_or_equal(0.0)
    end
  end

end
