facts("Generating set search") do

  ss = symmetric_search_space(3, (0.0, 1.0))
  @fact BlackBoxOptim.calc_initial_step_size(ss) --> (0.5 * (1.0 - 0.0))

  ss = symmetric_search_space(3, (-1.2, 42.0))
  @fact BlackBoxOptim.calc_initial_step_size(ss, 0.80) --> (0.80 * (42.0 + 1.2))  

end
