include("common.jl")

facts("Optimize single objective problems in 5, 10, 30 and 100 dimensions with sNES") do
  simple_problems = ["Sphere", "Schwefel2.22", "Schwefel2.22"]
  for problem in simple_problems
    context(problem) do
      p = BlackBoxOptim.Problems.examples[problem]

      @fact fitness_for_opt(p, 5, 20,  1e2, separable_nes) < 0.01 --> true
      @fact fitness_for_opt(p, 5, 20,  5e2, xnes) < 0.01 --> true

      @fact fitness_for_opt(p, 10, 20, 3e2, separable_nes) < 0.01 --> true
      @fact fitness_for_opt(p, 10, 20, 1e3, xnes) < 0.01 --> true

      @fact fitness_for_opt(p, 30, 25, 1e3, separable_nes) < 0.01 --> true
      @fact fitness_for_opt(p, 30, 25, 4e3, xnes) < 0.01 --> true

      @fact fitness_for_opt(p, 100, 25, 5e3, separable_nes) < 0.01 --> true
      # Cannot run xnes in 100 dimensions since it scales badly.
    end
  end

  context("Schwefel1.2") do
    problem = "Schwefel1.2"
    p = BlackBoxOptim.Problems.examples[problem]

    @fact fitness_for_opt(p, 30, 50, 4e3, separable_nes) < 10.0 --> true
    @fact fitness_for_opt(p, 30, 50, 8e3, xnes) < 10.0 --> true
  end

  context("Rosenbrock") do
    problem = "Rosenbrock"
    p = BlackBoxOptim.Problems.examples[problem]

    @fact fitness_for_opt(p, 30, 40, 2e4, separable_nes) < 100.0 --> true
    @fact fitness_for_opt(p, 30, 40, 5e4, xnes) < 100.0 --> true

    @fact fitness_for_opt(p, 50, 40, 3e4, separable_nes) < 100.0 --> true
  end

end
