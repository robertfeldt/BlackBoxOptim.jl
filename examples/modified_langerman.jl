# Code from Hans W. Borschers
function modLangerman(x)
    #if length(x) != 10; error("ModLangerman needs a 10-elements vector"); end
    sqr(z) = z^2
    x1 = x[1]; x2 = x[2]; x3 = x[3]; x4 = x[4]; x5 = x[5]
    x6 = x[6]; x7 = x[7]; x8 = x[8]; x9 = x[9]; x10 = x[10]
    
    y = 0.806*exp(-0.31830988618377*(sqr((-9.681) + x1) + sqr((-0.667) + x2) +
     sqr((-4.783) + x3) + sqr((-9.095) + x4) + sqr((-3.517) + x5) + sqr((-9.325
     ) + x6) + sqr((-6.544) + x7) + sqr((-0.211) + x8) + sqr((-5.122) + x9) +
     sqr((-2.02) + x10)))*cos(3.14159265359*(sqr((-9.681) + x1) + sqr((-0.667)
      + x2) + sqr((-4.783) + x3) + sqr((-9.095) + x4) + sqr((-3.517) + x5) +
     sqr((-9.325) + x6) + sqr((-6.544) + x7) + sqr((-0.211) + x8) + sqr((-5.122
     ) + x9) + sqr((-2.02) + x10))) + 0.517*exp(-0.31830988618377*(sqr((-9.4)
      + x1) + sqr((-2.041) + x2) + sqr((-3.788) + x3) + sqr((-7.931) + x4) +
     sqr((-2.882) + x5) + sqr((-2.672) + x6) + sqr((-3.568) + x7) + sqr((-1.284
     ) + x8) + sqr((-7.033) + x9) + sqr((-7.374) + x10)))*cos(3.14159265359*(
     sqr((-9.4) + x1) + sqr((-2.041) + x2) + sqr((-3.788) + x3) + sqr((-7.931)
      + x4) + sqr((-2.882) + x5) + sqr((-2.672) + x6) + sqr((-3.568) + x7) +
     sqr((-1.284) + x8) + sqr((-7.033) + x9) + sqr((-7.374) + x10))) + 0.1*exp(
     -0.31830988618377*(sqr((-8.025) + x1) + sqr((-9.152) + x2) + sqr((-5.114)
      + x3) + sqr((-7.621) + x4) + sqr((-4.564) + x5) + sqr((-4.711) + x6) +
     sqr((-2.996) + x7) + sqr((-6.126) + x8) + sqr((-0.734) + x9) + sqr((-4.982
     ) + x10)))*cos(3.14159265359*(sqr((-8.025) + x1) + sqr((-9.152) + x2) +
     sqr((-5.114) + x3) + sqr((-7.621) + x4) + sqr((-4.564) + x5) + sqr((-4.711
     ) + x6) + sqr((-2.996) + x7) + sqr((-6.126) + x8) + sqr((-0.734) + x9) +
     sqr((-4.982) + x10))) + 0.908*exp(-0.31830988618377*(sqr((-2.196) + x1) +
     sqr((-0.415) + x2) + sqr((-5.649) + x3) + sqr((-6.979) + x4) + sqr((-9.51)
      + x5) + sqr((-9.166) + x6) + sqr((-6.304) + x7) + sqr((-6.054) + x8) +
     sqr((-9.377) + x9) + sqr((-1.426) + x10)))*cos(3.14159265359*(sqr((-2.196)
      + x1) + sqr((-0.415) + x2) + sqr((-5.649) + x3) + sqr((-6.979) + x4) +
     sqr((-9.51) + x5) + sqr((-9.166) + x6) + sqr((-6.304) + x7) + sqr((-6.054)
      + x8) + sqr((-9.377) + x9) + sqr((-1.426) + x10))) + 0.965*exp(-
     0.31830988618377*(sqr((-8.074) + x1) + sqr((-8.777) + x2) + sqr((-3.467)
      + x3) + sqr((-1.867) + x4) + sqr((-6.708) + x5) + sqr((-6.349) + x6) +
     sqr((-4.534) + x7) + sqr((-0.276) + x8) + sqr((-7.633) + x9) + sqr((-1.567
     ) + x10)))*cos(3.14159265359*(sqr((-8.074) + x1) + sqr((-8.777) + x2) +
     sqr((-3.467) + x3) + sqr((-1.867) + x4) + sqr((-6.708) + x5) + sqr((-6.349
     ) + x6) + sqr((-4.534) + x7) + sqr((-0.276) + x8) + sqr((-7.633) + x9) +
     sqr((-1.567) + x10)))

    return -y
end

using BlackBoxOptim

mod_lang_problem = BlackBoxOptim.MinimizationProblemFamily(modLangerman, "ModifiedLangerman", 
  (0.0, 10.0), -0.9650)

BlackBoxOptim.repeated_bboptimize(5, mod_lang_problem, 10, [
  :generating_set_search, # GSS not very good on this one
  :adaptive_de_rand_1_bin_radiuslimited # While DE is at least better
  ],
  10.0, 1e-5, Dict{Symbol,Any}(
    :MinDeltaFitnessTolerance => 1e-50))
