using BlackBoxOptim

p = BlackBoxOptim.example_problems["Griewank"]

BlackBoxOptim.repeated_bboptimize(2, p, 30, [
  :generating_set_search, 
  #:random_search,
  #:separable_nes,
  :adaptive_de_rand_1_bin_radiuslimited],
  30.0, 1e-6, {:SaveFitnessTraceToCsv => true})