# Sphere function as stated in table II of the JADE paper:
#  http://150.214.190.154/EAMHCO/pdf/JADE.pdf
# Initial range: [-100, 100]^d
function sphere(x)
  sum(x.^2)
end

examples["Sphere"] = OptimizationProblem("Sphere",
                        sphere,
                        [-100, 100])