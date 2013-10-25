module SingleObjectiveProblems

  export OptimizationProblem, examples

  immutable OptimizationProblem
    name::ASCIIString
    f::Function
    initial_range::Vector{Float64}      
  end

  examples = Dict{ASCIIString, OptimizationProblem}()

  # Sphere function as stated in table II of the JADE paper:
  #  http://150.214.190.154/EAMHCO/pdf/JADE.pdf
  # Initial range: [-100, 100]^d
  function sphere(x)
    sum(x.^2)
  end

  examples["Sphere"] = OptimizationProblem("Sphere",
                        sphere,
                        [-100, 100])

  # Schwefel 2.22 function as stated in table II of the JADE paper:
  #  http://150.214.190.154/EAMHCO/pdf/JADE.pdf
  function schwefel2_22(x)
    ax = abs(x)
    sum(ax) + prod(ax)
  end
  
  examples["Schwefel2.22"] = OptimizationProblem("Schwefel2.22",
                              schwefel2_22,
                              [-10, 10])
  
  # Schwefel 1.2 function as stated in table II of the JADE paper:
  #  http://150.214.190.154/EAMHCO/pdf/JADE.pdf
  function schwefel1_2(x)
    s = 0
    for i in 1:length(x)
      partsum = sum(x[1:i])
      s+= (partsum^2)
    end
    s
  end
  
  examples["Schwefel1.2"] = OptimizationProblem("Schwefel1.2",
                              schwefel1_2,
                              [-100, 100])
  
  # Schwefel 2.21 function as stated in table II of the JADE paper:
  #  http://150.214.190.154/EAMHCO/pdf/JADE.pdf
  # Initial range: [-100, 100]^d
  function schwefel2_21(x)
    maximum(abs(x))
  end
  
  examples["Schwefel2.21"] = OptimizationProblem("Schwefel2.21",
                              schwefel2_22,
                              [-100, 100])

end # module SingleObjectiveProblems