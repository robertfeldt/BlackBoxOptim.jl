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

  context("hat_compare fitnesses of size 1 in a minimizing FitnessScheme") do
    scheme = float_vector_scheme_min()

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
    scheme = float_vector_scheme_max()

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
    scheme = float_vector_scheme_min()

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
    scheme = float_vector_scheme_max()

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

  context("isbetter/isworse/samefitness in a minimizing FitnessScheme") do
    scheme = float_vector_scheme_min()

    @fact isbetter([-1.0, 0.0], [1.0, 0.0], scheme) => true
    @fact isbetter([0.0, 0.0], [1.0, 0.0], scheme) => true

    @fact isbetter([1.0, 0.0], [-1.0, 0.0], scheme) => false
    @fact isbetter([1.0, 0.0], [0.0, 0.0], scheme) => false

    @fact isworse([-1.0, 0.0], [1.0, 0.0], scheme) => false
    @fact isworse([0.0, 0.0], [1.0, 0.0], scheme) => false

    @fact isworse([1.0, 0.0], [-1.0, 0.0], scheme) => true
    @fact isworse([1.0, 0.0], [0.0, 0.0], scheme) => true

    @fact samefitness([-1.0, 0.0], [1.0, 0.0], scheme) => false
    @fact samefitness([0.0, 0.0], [1.0, 0.0], scheme) => false
    @fact samefitness([0.0], [1.0], scheme) => false

    @fact samefitness([-1.0, 0.0], [-1.0, 0.0], scheme) => true
    @fact samefitness([0.0, 0.0], [0.0, 0.0], scheme) => true
    @fact samefitness([0.0], [0.0], scheme) => true
  end

  context("isbetter/isworse/samefitness in a maximizing FitnessScheme") do
    scheme = float_vector_scheme_max()

    @fact isbetter([-1.0, 0.0], [1.0, 0.0], scheme) => false
    @fact isbetter([0.0, 0.0], [1.0, 0.0], scheme) => false

    @fact isbetter([1.0, 0.0], [-1.0, 0.0], scheme) => true
    @fact isbetter([1.0, 0.0], [0.0, 0.0], scheme) => true

    @fact isworse([-1.0, 0.0], [1.0, 0.0], scheme) => true
    @fact isworse([0.0, 0.0], [1.0, 0.0], scheme) => true

    @fact isworse([1.0, 0.0], [-1.0, 0.0], scheme) => false
    @fact isworse([1.0, 0.0], [0.0, 0.0], scheme) => false

    @fact samefitness([-1.0, 0.0], [1.0, 0.0], scheme) => false
    @fact samefitness([0.0, 0.0], [1.0, 0.0], scheme) => false
    @fact samefitness([0.0], [1.0], scheme) => false

    @fact samefitness([-1.0, 0.0], [-1.0, 0.0], scheme) => true
    @fact samefitness([0.0, 0.0], [0.0, 0.0], scheme) => true
    @fact samefitness([0.0], [0.0], scheme) => true
  end

  context("worst_fitness") do
    @fact float_vector_scheme_max().worst_fitness => -Inf
    @fact float_vector_scheme_min().worst_fitness => Inf
  end

end