facts("Fitness") do
  context("hat_compare() Float64") do
    @fact hat_compare(1.0, 2.0) --> -1
    @fact hat_compare(-1.0, 1.0) --> -1

    @fact hat_compare(2.0, 1.0) --> 1
    @fact hat_compare(-2.0, -3.0) --> 1

    @fact hat_compare(1.0, 1.0) --> 0
    @fact hat_compare(0.0, 0.0) --> 0
    @fact hat_compare(-1.0, -1.0) --> 0
  end

  context("ScalarFitnessScheme") do
      context("is_minimizing()") do
        mins = MinimizingFitnessScheme
        @fact is_minimizing(mins) --> true

        maxs = MaximizingFitnessScheme
        @fact is_minimizing(maxs) --> false
      end

      context("hat_compare(..., MinimizingFitnessScheme)") do
        scheme = MinimizingFitnessScheme

        @fact hat_compare(1.0, 2.0, scheme) --> -1
        @fact hat_compare(2.0, 1.0, scheme) --> 1
        @fact hat_compare(1.0, 1.0, scheme) --> 0
      end

      context("hat_compare(..., MaximizingFitnessScheme)") do
        scheme = MaximizingFitnessScheme

        @fact hat_compare(1.0, 2.0, scheme) --> 1
        @fact hat_compare(2.0, 1.0, scheme) --> -1
        @fact hat_compare(1.0, 1.0, scheme) --> 0
      end

      # FIXME enable tests once v0.5 issue #14919 is fixed
      context("fitness_scheme(x, y)") do
        mins = MinimizingFitnessScheme
        @pending mins(5.0, 3.0) --> false
        @pending mins(1.0, 2.0) --> true

        maxs = MaximizingFitnessScheme
        @pending maxs(3.0, 5.0) --> false
        @pending maxs(2.0, 1.0) --> true
      end
  end

  context("TupleFitnessScheme") do
      context("general interface") do
        context("ParetoFitnessScheme{1}") do
            scheme = ParetoFitnessScheme{1}()
            @fact numobjectives(scheme) --> 1
            @fact is_minimizing(scheme) --> true
            @fact isnafitness(nafitness(scheme), scheme) --> true
            @fact isequal(nafitness(scheme), (NaN,)) --> true
            @fact isnafitness((NaN,), scheme) --> true
            @fact isnafitness((1.0,), scheme) --> false
        end

        context("ParetoFitnessScheme{1}(is_minimizing=false)") do
            scheme = ParetoFitnessScheme{1}(is_minimizing=false)
            @fact numobjectives(scheme) --> 1
            @fact is_minimizing(scheme) --> false
            @fact isnafitness(nafitness(scheme), scheme) --> true
            @fact isequal(nafitness(scheme), (NaN,)) --> true
            @fact isnafitness((NaN,), scheme) --> true
            @fact isnafitness((1.0,), scheme) --> false
        end

        context("ParetoFitnessScheme{3}()") do
            scheme = ParetoFitnessScheme{3}()
            @fact isequal(nafitness(scheme), (NaN,NaN,NaN)) --> true
            @fact numobjectives(scheme) --> 3
            @fact is_minimizing(scheme) --> true
            @fact isnafitness(nafitness(scheme), scheme) --> true
            @fact isnafitness((NaN,NaN,NaN), scheme) --> true
            @fact isnafitness((1.0,NaN,2.0), scheme) --> true
            @fact isnafitness((1.0,2.0,2.0), scheme) --> false
        end
      end

      context("hat_compare(..., ParetoFitnessScheme{1}())") do
        scheme = ParetoFitnessScheme{1}()

        @fact_throws hat_compare((-1.0,), (1.0, 2.0), scheme) MethodError

        @fact hat_compare((-1.0,), (1.0,), scheme) --> -1
        @fact hat_compare((0.0,), (0.3,), scheme) --> -1
        @fact hat_compare((11.3,), (354.65,), scheme) --> -1

        @fact hat_compare((-1.0,), (-1.0,), scheme) --> 0
        @fact hat_compare((0.0,), (0.0,), scheme) --> 0
        @fact hat_compare((0.2,), (0.2,), scheme) --> 0

        @fact hat_compare((1.0,), (0.0,), scheme) --> 1
        @fact hat_compare((0.0,), (-0.4,), scheme) --> 1
        @fact hat_compare((-0.65,), (-34.2,), scheme) --> 1
      end

      context("hat_compare(..., ParetoFitnessScheme{1}(is_minimizing=false))") do
        scheme = ParetoFitnessScheme{1}(is_minimizing=false)

        @fact hat_compare((-1.0,), (1.0,), scheme) --> 1
        @fact hat_compare((0.0,), (0.3,), scheme) --> 1
        @fact hat_compare((11.3,), (354.65,), scheme) --> 1

        @fact hat_compare((-1.0,), (-1.0,), scheme) --> 0
        @fact hat_compare((0.0,), (0.0,), scheme) --> 0
        @fact hat_compare((0.2,), (0.2,), scheme) --> 0

        @fact hat_compare((1.0,), (0.0,), scheme) --> -1
        @fact hat_compare((0.0,), (-0.4,), scheme) --> -1
        @fact hat_compare((-0.65,), (-34.2,), scheme) --> -1
      end

      context("hat_compare(..., ParetoFitnessScheme{2}(is_minimizing=true))") do
        scheme = ParetoFitnessScheme{2}(is_minimizing=true)
        @fact is_minimizing(scheme) --> true

        @fact_throws hat_compare((-1.0,), (1.0, 2.0), scheme) MethodError
        @fact_throws hat_compare((2.0, -1.0, 3.0), (1.0, 2.0), scheme) MethodError

        @fact hat_compare((-1.0, 0.0), (1.0, 0.0), scheme) --> -1
        @fact hat_compare((0.0, -1.0), (0.3, 1.0), scheme) --> -1

        @fact hat_compare((-1.0, 2.0), (-1.0, 2.0), scheme) --> 0
        @fact hat_compare((0.0, 0.0), (0.0, 0.0), scheme) --> 0
        @fact hat_compare((0.2, -0.4), (0.2, -0.4), scheme) --> 0
        @fact hat_compare((0.4, -0.4), (0.2, -0.2), scheme) --> 0
        @fact hat_compare((-0.65, 1.0), (1.0, -34.2), scheme) --> 0

        @fact hat_compare((11.3, 100.0), (11.3, 10.0), scheme) --> 1
        @fact hat_compare((1.0, 0.0), (0.0, -1.0), scheme) --> 1
        @fact hat_compare((0.0, 34.0), (-0.4, 20.0), scheme) --> 1
      end

      context("hat_compare(..., ParetoFitnessScheme{2}(is_minimizing=false))") do
        scheme = ParetoFitnessScheme{2}(is_minimizing=false)
        @fact is_minimizing(scheme) --> false

        @fact hat_compare((-1.0, 0.0), (1.0, 0.0), scheme) --> 1
        @fact hat_compare((0.0, -1.0), (0.3, 1.0), scheme) --> 1

        @fact hat_compare((-1.0, 2.0), (-1.0, 2.0), scheme) --> 0
        @fact hat_compare((0.0, 0.0), (0.0, 0.0), scheme) --> 0
        @fact hat_compare((0.2, -0.4), (0.2, -0.4), scheme) --> 0
        @fact hat_compare((0.4, -0.4), (0.2, -0.2), scheme) --> 0
        @fact hat_compare((-0.65, 1.0), (1.0, -34.2), scheme) --> 0

        @fact hat_compare((11.3, 100.0), (11.3, 10.0), scheme) --> -1
        @fact hat_compare((1.0, 0.0), (0.0, -1.0), scheme) --> -1
        @fact hat_compare((0.0, 34.0), (-0.4, 20.0), scheme) --> -1
      end

      context("is_better/is_worse/same_fitness(..., ParetoFitnessScheme{2}(is_minimizing=true))") do
        scheme = ParetoFitnessScheme{2}()

        @fact is_better((-1.0, 0.0), (1.0, 0.0), scheme) --> true
        @fact is_better((0.0, 0.0), (1.0, 0.0), scheme) --> true
        @fact_throws is_better((0.0,), (1.0,), scheme) MethodError

        @fact is_better((1.0, 0.0), (-1.0, 0.0), scheme) --> false
        @fact is_better((1.0, 0.0), (0.0, 0.0), scheme) --> false

        @fact is_worse((-1.0, 0.0), (1.0, 0.0), scheme) --> false
        @fact is_worse((0.0, 0.0), (1.0, 0.0), scheme) --> false
        @fact_throws is_worse((0.0, 1.0), (1.0,), scheme) MethodError

        @fact is_worse((1.0, 0.0), (-1.0, 0.0), scheme) --> true
        @fact is_worse((1.0, 0.0), (0.0, 0.0), scheme) --> true

        @fact same_fitness((-1.0, 0.0), (1.0, 0.0), scheme) --> false
        @fact same_fitness((0.0, 0.0), (1.0, 0.0), scheme) --> false
        @fact_throws same_fitness((0.0,), (1.0,), scheme) MethodError

        @fact same_fitness((-1.0, 0.0), (-1.0, 0.0), scheme) --> true
        @fact same_fitness((0.0, 0.0), (0.0, 0.0), scheme) --> true
        @fact_throws same_fitness((0.0,), (0.0,), scheme) MethodError
      end

      context("is_better/is_worse/same_fitness(..., ParetoFitnessScheme{2}(is_minimizing=false))") do
        scheme = ParetoFitnessScheme{2}(is_minimizing=false)

        @fact is_better((-1.0, 0.0), (1.0, 0.0), scheme) --> false
        @fact is_better((0.0, 0.0), (1.0, 0.0), scheme) --> false

        @fact is_better((1.0, 0.0), (-1.0, 0.0), scheme) --> true
        @fact is_better((1.0, 0.0), (0.0, 0.0), scheme) --> true

        @fact is_worse((-1.0, 0.0), (1.0, 0.0), scheme) --> true
        @fact is_worse((0.0, 0.0), (1.0, 0.0), scheme) --> true

        @fact is_worse((1.0, 0.0), (-1.0, 0.0), scheme) --> false
        @fact is_worse((1.0, 0.0), (0.0, 0.0), scheme) --> false

        @fact same_fitness((-1.0, 0.0), (1.0, 0.0), scheme) --> false
        @fact same_fitness((0.0, 0.0), (1.0, 0.0), scheme) --> false

        @fact same_fitness((-1.0, 0.0), (-1.0, 0.0), scheme) --> true
        @fact same_fitness((0.0, 0.0), (0.0, 0.0), scheme) --> true
      end
  end
  context("ϵ-box FitnessScheme") do
      context("ϵ-index()") do
        for (u, ϵ, is_minimizing, u_ix, delta) in [
                (2.3, 1.0, true, 2, 0.3),
                (2.3, 0.1, true, 23, 0.0),
                (2.3, 1.0, false, 3, 0.7),
                (1.0, 1.0, true, 1.0, 0.0),
                (2.0, 2.0, false, 1.0, 0.0),
                (-1.1, 1.0, true, -2, 0.9),
                (-5.4, 1.0, false, -5, 0.4)
            ]
          res = BlackBoxOptim.ϵ_index(u, ϵ, Val{is_minimizing})
          @fact res[1] --> u_ix
          @fact res[2] --> roughly(delta, atol=1E-12)
        end
      end

      context("IndexedTupleFitness") do
          fit = (0.55, 0.3, -0.21)
          ifit1 = IndexedTupleFitness(fit, sum(fit), 0.1, Val{true})
          @fact ifit1.orig --> fit
          @fact ifit1.agg --> sum(fit)
          @fact ifit1.index --> (5, 3, -3)
          @fact ifit1.dist --> roughly(norm([0.05, 0.09])/0.1)

          ifit2 = IndexedTupleFitness(fit, sum(fit), 0.1, Val{false})
          @fact ifit2.orig --> fit
          @fact ifit2.agg --> sum(fit)
          @fact ifit2.index --> (6, 3, -2)
          @fact ifit2.dist --> roughly(norm([0.05, 0.01])/0.1)

          ifit3 = convert(IndexedTupleFitness, fit,
                          EpsBoxDominanceFitnessScheme{3}(0.1))
          @fact ifit3 --> ifit1
          ifit4 = convert(IndexedTupleFitness, fit,
                          EpsBoxDominanceFitnessScheme{3}(0.1, is_minimizing=false))
          @fact ifit4 --> ifit2
      end

      context("hat_compare(..., EpsBoxDominanceFitnessScheme{2}(...))") do
        minscheme = EpsBoxDominanceFitnessScheme{2}(0.1, is_minimizing=true)
        maxscheme = EpsBoxDominanceFitnessScheme{2}(0.1, is_minimizing=false)

        @fact isnafitness(nafitness(minscheme), minscheme) --> true
        @fact isnafitness(nafitness(maxscheme), minscheme) --> true
        @fact isnafitness(nafitness(minscheme), maxscheme) --> true
        @fact isnafitness(nafitness(maxscheme), maxscheme) --> true

        @fact hat_compare((-1.0, 0.0), (-1.0, 0.0), minscheme) --> (0, true)
        @fact hat_compare((-1.0, 0.0), (-1.0, 0.0), maxscheme) --> (0, true)
        @fact hat_compare((-1.0, 0.0), (1.0, 0.0), minscheme) --> (-1, false)
        @fact hat_compare((-1.0, 0.0), (1.0, 0.0), maxscheme) --> (1, false)
        @fact hat_compare((-1.0, 0.0), (-0.99, 0.0), minscheme) --> (-1, true)
        @fact hat_compare((-1.0, 0.0), (-0.99, 0.0), maxscheme) --> (1, false)
        @fact hat_compare((-1.0, 0.0), (-1.01, 0.0), minscheme) --> (1, false)
        @fact hat_compare((-1.0, 0.0), (-1.01, 0.0), maxscheme) --> (-1, true)

        @fact hat_compare((-0.97, 0.44), (-0.96, 0.43), minscheme) --> (0, true)
        @fact hat_compare((-0.97, 0.44), (-0.96, 0.43), maxscheme) --> (0, true)
        @fact hat_compare((-1.0, 0.5), (-1.0, 0.51), minscheme) --> (-1, true)
        @fact hat_compare((-1.0, 0.5), (-1.0, 0.51), maxscheme) --> (1, false)

        @fact hat_compare((0.9, 0.0), (1.0, 0.0), minscheme) --> (-1, false)
        @fact hat_compare((0.9, 0.0), (1.0, 0.0), maxscheme) --> (1, false)
        @fact hat_compare((0.95, 0.0), (0.98, 0.0), minscheme) --> (-1, true)
        @fact hat_compare((0.95, 0.0), (0.98, 0.0), maxscheme) --> (1, true)

        @fact hat_compare((-1.0, 0.05), (1.0, 0.0), minscheme) --> (-1, false)
        @fact hat_compare((-1.0, 0.05), (1.0, 0.0), maxscheme) --> (0, false)
        @fact hat_compare((0.9, 0.05), (1.0, 0.0), minscheme) --> (-1, false)
        @fact hat_compare((0.9, 0.05), (1.0, 0.0), maxscheme) --> (0, false)
        @fact hat_compare((0.95, 0.05), (0.98, 0.0), minscheme) --> (-1, true)
        @fact hat_compare((0.95, 0.05), (0.98, 0.0), maxscheme) --> (-1, false)
        @fact hat_compare((0.95, 0.05), (0.98, 0.01), maxscheme) --> (-1, true)
        @fact hat_compare((0.95, 0.05), (0.96, 0.0), minscheme) --> (1, true)
        @fact hat_compare((0.95, 0.05), (0.96, 0.0), maxscheme) --> (-1, false)
        @fact hat_compare((0.95, 0.05), (0.96, 0.04), minscheme) --> (-1, true)
        @fact hat_compare((0.95, 0.05), (0.96, 0.04), maxscheme) --> (-1, true)

        @fact hat_compare((6.9928266604943286,0.03386770153536918),
                          (7.007808609410211,0.032833634435236035), minscheme, 0) --> (-1, false)
        @fact hat_compare((6.9928266604943286,0.03386770153536918),
                          (7.007808609410211,0.032833634435236035), minscheme, 1) --> (-1, false)
        @fact hat_compare((6.9928266604943286,0.03386770153536918),
                          (7.007808609410211,0.032833634435236035), minscheme, -1) --> (-1, false)
      end
  end
end
