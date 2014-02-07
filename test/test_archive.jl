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
    @fact delta_fitness(a)    => Inf

    add_candidate!(a, 2.0, [1.0])
    @fact a.size            => 3
    @fact a.count           => 2
    @fact best_fitness(a)   => 1.0
    @fact best_candidate(a) => [0.0]
    @fact last_top_fitness(a) => 2.0
    @fact delta_fitness(a)    => Inf

    add_candidate!(a, 0.5, [2.0])
    @fact a.size            => 3
    @fact a.count           => 3
    @fact best_fitness(a)   => 0.5
    @fact best_candidate(a) => [2.0]
    @fact last_top_fitness(a) => 2.0
    @fact delta_fitness(a)    => 0.5

    add_candidate!(a, 0.8, [4.0])
    @fact a.size            => 3
    @fact a.count           => 3
    @fact best_fitness(a)   => 0.5
    @fact best_candidate(a) => [2.0]
    @fact last_top_fitness(a) => 1.0

    expected = ((0.8 - 0.5) / ((1 - 0.05)^(-2/1) - 1))
    @fact width_of_confidence_interval(a, 0.05) => expected
    @fact fitness_improvement_potential(a, 0.05) => (expected / 0.50)

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

  context("magnitude_class for positive fitness values") do

    @fact BlackBoxOptim.magnitude_class(1.0)      => (1.0, 0.0)
    @fact BlackBoxOptim.magnitude_class(10.0)     => (1.0, 1.0)
    @fact BlackBoxOptim.magnitude_class(1e2)      => (1.0, 2.0)
    @fact BlackBoxOptim.magnitude_class(1e102)    => (1.0, 102.0)
    @fact BlackBoxOptim.magnitude_class(1e-1)     => (1.0, -1.0)
    @fact BlackBoxOptim.magnitude_class(1e-2)     => (1.0, -2.0)
    @fact BlackBoxOptim.magnitude_class(1e-9)     => (1.0, -9.0)
    @fact BlackBoxOptim.magnitude_class(1e-12)    => (1.0, -12.0)

    @fact BlackBoxOptim.magnitude_class(2.0)      => (1.0, 0.3)
    @fact BlackBoxOptim.magnitude_class(20.0)     => (1.0, 1.3)
    @fact BlackBoxOptim.magnitude_class(200.0)    => (1.0, 2.3)

    @fact BlackBoxOptim.magnitude_class(5.0)      => (1.0, 0.6)
    @fact BlackBoxOptim.magnitude_class(50.0)     => (1.0, 1.6)
    @fact BlackBoxOptim.magnitude_class(500.0)    => (1.0, 2.6)

  end

  context("magnitude_class for negative fitness values") do

    @fact BlackBoxOptim.magnitude_class(-1.0)      => (-1.0, 0.0)
    @fact BlackBoxOptim.magnitude_class(-10.0)     => (-1.0, 1.0)
    @fact BlackBoxOptim.magnitude_class(-1e2)      => (-1.0, 2.0)
    @fact BlackBoxOptim.magnitude_class(-1e102)    => (-1.0, 102.0)
    @fact BlackBoxOptim.magnitude_class(-1e-1)     => (-1.0, -1.0)
    @fact BlackBoxOptim.magnitude_class(-1e-2)     => (-1.0, -2.0)
    @fact BlackBoxOptim.magnitude_class(-1e-9)     => (-1.0, -9.0)
    @fact BlackBoxOptim.magnitude_class(-1e-12)    => (-1.0, -12.0)

    @fact BlackBoxOptim.magnitude_class(-2.0)      => (-1.0, 0.3)
    @fact BlackBoxOptim.magnitude_class(-20.0)     => (-1.0, 1.3)
    @fact BlackBoxOptim.magnitude_class(-200.0)    => (-1.0, 2.3)

    @fact BlackBoxOptim.magnitude_class(-5.0)      => (-1.0, 0.6)
    @fact BlackBoxOptim.magnitude_class(-50.0)     => (-1.0, 1.6)
    @fact BlackBoxOptim.magnitude_class(-500.0)    => (-1.0, 2.6)

  end

end