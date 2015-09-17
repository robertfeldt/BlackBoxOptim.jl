facts("TopListArchive") do
  context("ArchivedIndividual") do
    # test equality
    @fact isequal(BlackBoxOptim.ArchivedIndividual([3.0, 2.0], 2.0), BlackBoxOptim.ArchivedIndividual([3.0, 2.0], 2.0)) --> true
    @fact (BlackBoxOptim.ArchivedIndividual([3.0, 2.0], 2.0) == BlackBoxOptim.ArchivedIndividual([3.0, 2.0], 2.0)) --> true
    @fact (BlackBoxOptim.ArchivedIndividual([3.0, 2.0], 2.0) != BlackBoxOptim.ArchivedIndividual([3.0, 2.0], 2.0)) --> false
    @fact (BlackBoxOptim.ArchivedIndividual([3.0, 2.0], 1.0) != BlackBoxOptim.ArchivedIndividual([3.0, 2.0], 2.0)) --> true
    @fact (BlackBoxOptim.ArchivedIndividual([1.0, 2.0], 2.0) != BlackBoxOptim.ArchivedIndividual([3.0, 2.0], 2.0)) --> true
  end

  context("Constructing a small archive and adding to it") do

    a = TopListArchive(MinimizingFitnessScheme, 1, 3)

    @fact capacity(a) --> 3
    @fact length(a)   --> 0
    @fact best_fitness(a) --> isnan

    BlackBoxOptim.add_candidate!(a, 1.0, [0.0])
    @fact capacity(a)         --> 3
    @fact length(a)           --> 1
    @fact best_fitness(a)     --> 1.0
    @fact best_candidate(a)   --> [0.0]
    @fact last_top_fitness(a) --> 1.0
    @fact delta_fitness(a)    --> Inf

    BlackBoxOptim.add_candidate!(a, 2.0, [1.0])
    @fact capacity(a)       --> 3
    @fact length(a)         --> 2
    @fact best_fitness(a)   --> 1.0
    @fact best_candidate(a) --> [0.0]
    @fact last_top_fitness(a) --> 2.0
    @fact delta_fitness(a)    --> Inf

    BlackBoxOptim.add_candidate!(a, 0.5, [2.0])
    @fact capacity(a)         --> 3
    @fact length(a)           --> 3
    @fact best_fitness(a)   --> 0.5
    @fact best_candidate(a) --> [2.0]
    @fact last_top_fitness(a) --> 2.0
    @fact delta_fitness(a)    --> 0.5

    BlackBoxOptim.add_candidate!(a, 0.8, [4.0])
    @fact capacity(a)         --> 3
    @fact length(a)           --> 3
    @fact best_fitness(a)   --> 0.5
    @fact best_candidate(a) --> [2.0]
    @fact last_top_fitness(a) --> 1.0
    @fact delta_fitness(a)    --> 0.5

    expected = ((0.8 - 0.5) / ((1 - 0.05)^(-2/1) - 1))
    @fact width_of_confidence_interval(a, 0.05) --> expected
    @fact fitness_improvement_potential(a, 0.05) --> (expected / 0.50)

    BlackBoxOptim.add_candidate!(a, 0.4, [1.9])
    @fact capacity(a)       --> 3
    @fact length(a)         --> 3
    @fact best_fitness(a)   --> 0.4
    @fact best_candidate(a) --> [1.9]
    @fact last_top_fitness(a) --> 0.8
    @fact delta_fitness(a)    --> roughly(0.1)

    expected = ((0.5 - 0.4) / ((1 - 0.05)^(-2/1) - 1))
    @fact width_of_confidence_interval(a, 0.05) --> expected
    @fact fitness_improvement_potential(a, 0.05) --> (expected / 0.40)

    # identical candidate is not inserted
    BlackBoxOptim.add_candidate!(a, 0.5, [2.0])
    @fact capacity(a)       --> 3
    @fact length(a)         --> 3
    @fact last_top_fitness(a)   --> 0.8

    # different candidates with the same fitness are inserted
    BlackBoxOptim.add_candidate!(a, 0.5, [3.0])
    @fact length(a)             --> 3
    @fact last_top_fitness(a)   --> 0.5
  end

  context("magnitude_class for positive fitness values") do

    @fact BlackBoxOptim.magnitude_class(1.0)      --> (1.0, 0.0)
    @fact BlackBoxOptim.magnitude_class(10.0)     --> (1.0, 1.0)
    @fact BlackBoxOptim.magnitude_class(1e2)      --> (1.0, 2.0)
    @fact BlackBoxOptim.magnitude_class(1e102)    --> (1.0, 102.0)
    @fact BlackBoxOptim.magnitude_class(1e-1)     --> (1.0, -1.0)
    @fact BlackBoxOptim.magnitude_class(1e-2)     --> (1.0, -2.0)
    @fact BlackBoxOptim.magnitude_class(1e-9)     --> (1.0, -9.0)
    @fact BlackBoxOptim.magnitude_class(1e-12)    --> (1.0, -12.0)

    @fact BlackBoxOptim.magnitude_class(2.0)      --> (1.0, 0.3)
    @fact BlackBoxOptim.magnitude_class(20.0)     --> (1.0, 1.3)
    @fact BlackBoxOptim.magnitude_class(200.0)    --> (1.0, 2.3)

    @fact BlackBoxOptim.magnitude_class(5.0)      --> (1.0, 0.6)
    @fact BlackBoxOptim.magnitude_class(50.0)     --> (1.0, 1.6)
    @fact BlackBoxOptim.magnitude_class(500.0)    --> (1.0, 2.6)

  end

  context("magnitude_class for negative fitness values") do

    @fact BlackBoxOptim.magnitude_class(-1.0)      --> (-1.0, 0.0)
    @fact BlackBoxOptim.magnitude_class(-10.0)     --> (-1.0, 1.0)
    @fact BlackBoxOptim.magnitude_class(-1e2)      --> (-1.0, 2.0)
    @fact BlackBoxOptim.magnitude_class(-1e102)    --> (-1.0, 102.0)
    @fact BlackBoxOptim.magnitude_class(-1e-1)     --> (-1.0, -1.0)
    @fact BlackBoxOptim.magnitude_class(-1e-2)     --> (-1.0, -2.0)
    @fact BlackBoxOptim.magnitude_class(-1e-9)     --> (-1.0, -9.0)
    @fact BlackBoxOptim.magnitude_class(-1e-12)    --> (-1.0, -12.0)

    @fact BlackBoxOptim.magnitude_class(-2.0)      --> (-1.0, 0.3)
    @fact BlackBoxOptim.magnitude_class(-20.0)     --> (-1.0, 1.3)
    @fact BlackBoxOptim.magnitude_class(-200.0)    --> (-1.0, 2.3)

    @fact BlackBoxOptim.magnitude_class(-5.0)      --> (-1.0, 0.6)
    @fact BlackBoxOptim.magnitude_class(-50.0)     --> (-1.0, 1.6)
    @fact BlackBoxOptim.magnitude_class(-500.0)    --> (-1.0, 2.6)

  end

  context("archive copies the individuals") do
    a = TopListArchive(MinimizingFitnessScheme, 2, 3)

    indiv = [0.0, 2.0]
    BlackBoxOptim.add_candidate!(a, 1.0, indiv)
    # equal but not identical
    @fact (best_candidate(a) == indiv) --> true
    @fact (best_candidate(a) === indiv) --> false

    # modify the vector
    indiv[1] = 5.0
    # test that archived version is still the same
    @fact best_candidate(a) --> [0.0, 2.0]
  end

  context("for maximizing fitness") do
    a = TopListArchive(MaximizingFitnessScheme, 1, 3)

    BlackBoxOptim.add_candidate!(a, 1.0, [0.0])
    @fact best_fitness(a)     --> 1.0
    @fact best_candidate(a)   --> [0.0]
    @fact last_top_fitness(a) --> 1.0
    @fact delta_fitness(a)    --> Inf

    BlackBoxOptim.add_candidate!(a, 2.0, [1.0])
    @fact best_fitness(a)   --> 2.0
    @fact best_candidate(a) --> [1.0]
    @fact last_top_fitness(a) --> 1.0
    @fact delta_fitness(a)    --> 1.0

    BlackBoxOptim.add_candidate!(a, 0.5, [2.0])
    @fact best_fitness(a)   --> 2.0
    @fact best_candidate(a) --> [1.0]
    @fact last_top_fitness(a) --> 0.5
    @fact delta_fitness(a)    --> 1.0

    BlackBoxOptim.add_candidate!(a, 1.5, [4.0])
    @fact best_fitness(a)   --> 2.0
    @fact best_candidate(a) --> [1.0]
    @fact last_top_fitness(a) --> 1.0
    @fact delta_fitness(a)    --> 1.0

    BlackBoxOptim.add_candidate!(a, 2.5, [5.0])
    @fact best_fitness(a)   --> 2.5
    @fact best_candidate(a) --> [5.0]
    @fact last_top_fitness(a) --> 1.5
    @fact delta_fitness(a)    --> 0.5
  end
end
