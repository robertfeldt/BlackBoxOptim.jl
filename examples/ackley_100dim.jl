using BlackBoxOptim

p = BlackBoxOptim.example_problems["Ackley"]

BlackBoxOptim.repeated_bboptimize(5, p, 100, [
  :generating_set_search, 
  :adaptive_de_rand_1_bin_radiuslimited,
  :random_search,
  :separable_nes], 
  20.0)