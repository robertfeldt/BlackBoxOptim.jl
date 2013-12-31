facts("TopListArchive") do

  context("Constructing a small archive and adding to it") do

    a = TopListArchive(1, 3)

    @fact a.size  => 3
    @fact a.count => 0

    add_candidate!(a, 1.0, [0.0])
    @fact a.size              => 3
    @fact a.count             => 1
    @fact best_fitness(a)     => 1.0
    @fact best_candidate(a)   => [0.0]
    @fact last_top_fitness(a) => 1.0

    add_candidate!(a, 2.0, [1.0])
    @fact a.size            => 3
    @fact a.count           => 2
    @fact best_fitness(a)   => 1.0
    @fact best_candidate(a) => [0.0]
    @fact last_top_fitness(a) => 2.0

    add_candidate!(a, 0.5, [2.0])
    @fact a.size            => 3
    @fact a.count           => 3
    @fact best_fitness(a)   => 0.5
    @fact best_candidate(a) => [2.0]
    @fact last_top_fitness(a) => 2.0

    add_candidate!(a, 0.8, [4.0])
    @fact a.size            => 3
    @fact a.count           => 3
    @fact best_fitness(a)   => 0.5
    @fact best_candidate(a) => [2.0]
    @fact last_top_fitness(a) => 1.0

    add_candidate!(a, 0.4, [1.9])
    @fact a.size            => 3
    @fact a.count           => 3
    @fact best_fitness(a)   => 0.4
    @fact best_candidate(a) => [1.9]
    @fact last_top_fitness(a) => 0.8

    expected = ((0.5 - 0.4) / ((1 - 0.05)^(-2/1) - 1))
    @fact width_of_confidence_interval(a, 0.05) => expected
    @fact fitness_improvement_potential(a, 0.05) => (expected / 0.40)
  end

end