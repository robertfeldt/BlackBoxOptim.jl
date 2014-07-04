function rosenbrock2d(x)
  return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
end

facts("bboptimize smoketest") do

 context("test each method option in short runs, just to ensure they work") do

   for(m in keys(BlackBoxOptim.ValidMethods))
     b, f = bboptimize(rosenbrock2d; method = m, 
      search_space = [(-5.0, 5.0), (-2.0, 2.0)], max_time = 0.3,
      parameters = {:ShowTrace => false})
     if (m != :random_search) & (m != :de_rand_2_bin)
       @fact f < 10.0 => true
     end
   end

 end

end