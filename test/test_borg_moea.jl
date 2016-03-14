function schaffer1(x)
  return sumabs2(x), sumabs2(x .- 2.0)
end

facts("BorgMOEA") do
  res = bboptimize(schaffer1; Method=:borg_moea,
                   FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true),
                   SearchRange=(-10.0, 10.0), NumDimensions=2, ϵ=0.01,
                   MaxSteps=5000, TraceMode=:silent)
end
