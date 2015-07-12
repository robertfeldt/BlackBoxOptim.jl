facts("Fitness") do

  context("hat_compare floats") do
    @fact hat_compare(1.0, 2.0) => -1
    @fact hat_compare(-1.0, 1.0) => -1

    @fact hat_compare(2.0, 1.0) => 1
    @fact hat_compare(-2.0, -3.0) => 1

    @fact hat_compare(1.0, 1.0) => 0
    @fact hat_compare(0.0, 0.0) => 0
    @fact hat_compare(-1.0, -1.0) => 0
  end

  context("is_minimizing in ScalarFitness schemes") do
    mins = BlackBoxOptim.ScalarFitness{true}()
    @fact BlackBoxOptim.is_minimizing(mins) => true

    maxs = BlackBoxOptim.ScalarFitness{false}()
    @fact BlackBoxOptim.is_minimizing(maxs) => false
  end

  context("hat_compare floats in a minimizing fitness scheme") do
    scheme = BlackBoxOptim.ScalarFitness{true}()

    @fact hat_compare(1.0, 2.0, scheme) => -1
    @fact hat_compare(2.0, 1.0, scheme) => 1
    @fact hat_compare(1.0, 1.0, scheme) => 0
  end

  context("hat_compare floats in a maximizing fitness scheme") do
    scheme = BlackBoxOptim.ScalarFitness{false}()

    @fact hat_compare(1.0, 2.0, scheme) => 1
    @fact hat_compare(2.0, 1.0, scheme) => -1
    @fact hat_compare(1.0, 1.0, scheme) => 0
  end

  context("hat_compare fitnesses of size 1 in a minimizing FitnessScheme") do
    scheme = vector_fitness_scheme_min(1)
    @fact BlackBoxOptim.is_minimizing(scheme) => true

    @fact hat_compare([-1.0], [1.0], scheme) => -1
    @fact hat_compare([0.0], [0.3], scheme) => -1
    @fact hat_compare([11.3], [354.65], scheme) => -1

    @fact hat_compare([-1.0], [-1.0], scheme) => 0
    @fact hat_compare([0.0], [0.0], scheme) => 0
    @fact hat_compare([0.2], [0.2], scheme) => 0

    @fact hat_compare([1.0], [0.0], scheme) => 1
    @fact hat_compare([0.0], [-0.4], scheme) => 1
    @fact hat_compare([-0.65], [-34.2], scheme) => 1
  end

  context("hat_compare fitnesses of size 1 in a maximizing FitnessScheme") do
    scheme = vector_fitness_scheme_max(1)
    @fact BlackBoxOptim.is_minimizing(scheme) => false

    @fact hat_compare([-1.0], [1.0], scheme) => 1
    @fact hat_compare([0.0], [0.3], scheme) => 1
    @fact hat_compare([11.3], [354.65], scheme) => 1

    @fact hat_compare([-1.0], [-1.0], scheme) => 0
    @fact hat_compare([0.0], [0.0], scheme) => 0
    @fact hat_compare([0.2], [0.2], scheme) => 0

    @fact hat_compare([1.0], [0.0], scheme) => -1
    @fact hat_compare([0.0], [-0.4], scheme) => -1
    @fact hat_compare([-0.65], [-34.2], scheme) => -1
  end

  context("hat_compare fitnesses of size > 1 in a minimizing FitnessScheme") do
    scheme = vector_fitness_scheme_min(2)
    @fact BlackBoxOptim.is_minimizing(scheme) => true

    @fact hat_compare([-1.0, 0.0], [1.0, 0.0], scheme) => -1
    @fact hat_compare([0.0, -1.0], [0.3, 1.0], scheme) => -1
    @fact hat_compare([11.3, 100.0], [354.65, 10.0], scheme) => -1

    @fact hat_compare([-1.0, 2.0], [-1.0, 2.0], scheme) => 0
    @fact hat_compare([0.0, 0.0], [0.0, 0.0], scheme) => 0
    @fact hat_compare([0.2, -0.4], [0.2, -0.4], scheme) => 0

    @fact hat_compare([1.0, 0.0], [0.0, -1.0], scheme) => 1
    @fact hat_compare([0.0, 34.0], [-0.4, 20.0], scheme) => 1
    @fact hat_compare([-0.65, 1.0], [1.0, -34.2], scheme) => 1
  end

  context("hat_compare fitnesses of size > 1 in a minimizing FitnessScheme") do
    scheme = vector_fitness_scheme_max(2)
    @fact BlackBoxOptim.is_minimizing(scheme) => false

    @fact hat_compare([-1.0, 0.0], [1.0, 0.0], scheme) => 1
    @fact hat_compare([0.0, -1.0], [0.3, 1.0], scheme) => 1
    @fact hat_compare([11.3, 100.0], [354.65, 10.0], scheme) => 1

    @fact hat_compare([-1.0, 2.0], [-1.0, 2.0], scheme) => 0
    @fact hat_compare([0.0, 0.0], [0.0, 0.0], scheme) => 0
    @fact hat_compare([0.2, -0.4], [0.2, -0.4], scheme) => 0

    @fact hat_compare([1.0, 0.0], [0.0, -1.0], scheme) => -1
    @fact hat_compare([0.0, 34.0], [-0.4, 20.0], scheme) => -1
    @fact hat_compare([-0.65, 1.0], [1.0, -34.2], scheme) => -1
  end

  context("is_better/is_worse/same_fitness in a minimizing FitnessScheme") do
    scheme = vector_fitness_scheme_min(2)

    @fact is_better([-1.0, 0.0], [1.0, 0.0], scheme) => true
    @fact is_better([0.0, 0.0], [1.0, 0.0], scheme) => true

    @fact is_better([1.0, 0.0], [-1.0, 0.0], scheme) => false
    @fact is_better([1.0, 0.0], [0.0, 0.0], scheme) => false

    @fact is_worse([-1.0, 0.0], [1.0, 0.0], scheme) => false
    @fact is_worse([0.0, 0.0], [1.0, 0.0], scheme) => false

    @fact is_worse([1.0, 0.0], [-1.0, 0.0], scheme) => true
    @fact is_worse([1.0, 0.0], [0.0, 0.0], scheme) => true

    @fact same_fitness([-1.0, 0.0], [1.0, 0.0], scheme) => false
    @fact same_fitness([0.0, 0.0], [1.0, 0.0], scheme) => false
    @fact same_fitness([0.0], [1.0], scheme) => false

    @fact same_fitness([-1.0, 0.0], [-1.0, 0.0], scheme) => true
    @fact same_fitness([0.0, 0.0], [0.0, 0.0], scheme) => true
    @fact same_fitness([0.0], [0.0], scheme) => true
  end

  context("is_better/is_worse/same_fitness in a maximizing FitnessScheme") do
    scheme = vector_fitness_scheme_max(2)

    @fact is_better([-1.0, 0.0], [1.0, 0.0], scheme) => false
    @fact is_better([0.0, 0.0], [1.0, 0.0], scheme) => false

    @fact is_better([1.0, 0.0], [-1.0, 0.0], scheme) => true
    @fact is_better([1.0, 0.0], [0.0, 0.0], scheme) => true

    @fact is_worse([-1.0, 0.0], [1.0, 0.0], scheme) => true
    @fact is_worse([0.0, 0.0], [1.0, 0.0], scheme) => true

    @fact is_worse([1.0, 0.0], [-1.0, 0.0], scheme) => false
    @fact is_worse([1.0, 0.0], [0.0, 0.0], scheme) => false

    @fact same_fitness([-1.0, 0.0], [1.0, 0.0], scheme) => false
    @fact same_fitness([0.0, 0.0], [1.0, 0.0], scheme) => false
    @fact same_fitness([0.0], [1.0], scheme) => false

    @fact same_fitness([-1.0, 0.0], [-1.0, 0.0], scheme) => true
    @fact same_fitness([0.0, 0.0], [0.0, 0.0], scheme) => true
    @fact same_fitness([0.0], [0.0], scheme) => true
  end

  if VERSION >= v"0.4.0-dev+1258" # FIXME remove version check once v0.4 is released
    context("fitness_scheme(x, y)") do
      mins = BlackBoxOptim.ScalarFitness{true}()
      @fact mins(5.0, 3.0) => false
      @fact mins(1.0, 2.0) => true

      maxs = BlackBoxOptim.ScalarFitness{false}()
      @fact maxs(3.0, 5.0) => false
      @fact maxs(2.0, 1.0) => true
    end
  end
end
