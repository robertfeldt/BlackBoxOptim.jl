using BlackBoxOptim

problem = BlackBoxOptim.example_problems["Ackley"]
numrepeats = 5
dim = 100
methods = [
  :generating_set_search, 
  :probabilistic_descent, 
  :adaptive_de_rand_1_bin_radiuslimited,
  :random_search,
  :separable_nes]
max_time = 10.0
ftol = 1e-5

BlackBoxOptim.repeated_bboptimize(numrepeats, problem, dim, methods, max_time, ftol)