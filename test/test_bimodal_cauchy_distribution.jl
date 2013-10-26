facts("Bimodal Cauchy Distributions") do

  context("sample bimodal cauchy with truncation on one side") do
    bc = GlobalOptim.bimodal_cauchy(0.65, 0.1, 1.0, 0.1)

    for(i in 1:100)
      v = GlobalOptim.sample_bimodal_cauchy(bc; truncateAbove1 = true, truncateBelow0 = false)

      @fact v <= 1.0 => true
      @fact v > 0.0 => true # Very unlikely to be 0.0!!?
    end
  end

  context("sample bimodal cauchy with truncation on both sides") do
    bc = GlobalOptim.bimodal_cauchy(0.1, 0.1, 0.95, 0.1)

    for(i in 1:100)
      v = GlobalOptim.sample_bimodal_cauchy(bc; truncateAbove1 = true, truncateBelow0 = true)

      @fact v <= 1.0 => true
      @fact v >= 0.0 => true
    end
  end

end