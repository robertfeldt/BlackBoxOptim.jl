facts("EpsBoxArchive") do
  scheme = EpsBoxDominanceFitnessScheme{2}(0.1, is_minimizing=true)

  context("EpsBoxFrontierIndividual") do
    arch_indiv = BlackBoxOptim.EpsBoxFrontierIndividual(
            convert(IndexedTupleFitness, (2.21, 1.12), scheme), [3.0, 2.0], 3, 5, 2)
    @fact arch_indiv.fitness.orig --> (2.21, 1.12)
    @fact arch_indiv.fitness.index --> (22, 11)
    @fact arch_indiv.timestamp --> greater_than(0)
    @fact arch_indiv.tag --> 3
    @fact arch_indiv.num_fevals --> 5
    @fact arch_indiv.n_restarts --> 2
  end

  context("Constructing a small archive and adding to it") do
    a = EpsBoxArchive(scheme, max_size=100)

    @fact capacity(a) --> 100
    @fact length(a)   --> 0
    @fact BlackBoxOptim.candidates_without_progress(a) --> 0
    @fact BlackBoxOptim.tagcounts(a) --> Dict{Int,Int}()

    BlackBoxOptim.add_candidate!(a, convert(IndexedTupleFitness, (2.21, 1.12), scheme), [0.0, 5.0], 1, 8)
    @fact capacity(a)         --> 100
    @fact length(a)           --> 1
    @fact a.best_candidate_ix --> 1
    @fact a.frontier[1].num_fevals --> 8
    @fact BlackBoxOptim.candidates_without_progress(a) --> 0
    @fact BlackBoxOptim.tagcounts(a) --> Dict{Int,Int}(1=>1)

    BlackBoxOptim.add_candidate!(a, (3.21, 1.12), [1.0, 5.0], 2)
    @fact capacity(a)         --> 100
    @fact length(a)           --> 1
    @fact a.best_candidate_ix --> 1
    @fact a.frontier[1].num_fevals --> 8
    @fact BlackBoxOptim.candidates_without_progress(a) --> 1
    @fact BlackBoxOptim.tagcounts(a) --> Dict{Int,Int}(1=>1)

    BlackBoxOptim.add_candidate!(a, (1.21, 2.11), [1.2, 3.0], 3)
    @fact capacity(a)         --> 100
    @fact length(a)           --> 2
    @fact a.frontier[2].num_fevals --> 3
    @fact a.best_candidate_ix --> 2
    @fact BlackBoxOptim.candidates_without_progress(a) --> 0
    @fact BlackBoxOptim.tagcounts(a) --> Dict{Int,Int}(1=>1, 3=>1)

    BlackBoxOptim.add_candidate!(a, (0.52, 3.15), [1.5, 3.0], 3)
    @fact capacity(a)         --> 100
    @fact length(a)           --> 3
    @fact a.best_candidate_ix --> 2
    @fact BlackBoxOptim.candidates_without_progress(a) --> 0
    @fact BlackBoxOptim.tagcounts(a) --> Dict{Int,Int}(1=>1, 3=>2)

    BlackBoxOptim.add_candidate!(a, (0.21, 0.21), [1.5, 3.0], 5)
    @fact capacity(a)         --> 100
    @fact length(a)           --> 1
    @fact a.best_candidate_ix --> 1
    @fact BlackBoxOptim.candidates_without_progress(a) --> 0
    @fact BlackBoxOptim.tagcounts(a) --> Dict{Int,Int}(5=>1)

    BlackBoxOptim.add_candidate!(a, (0.2, 0.2), [1.5, 3.0], 6)
    @fact capacity(a)         --> 100
    @fact length(a)           --> 1
    @fact a.best_candidate_ix --> 1
    @fact BlackBoxOptim.candidates_without_progress(a) --> 1
    @fact BlackBoxOptim.tagcounts(a) --> Dict{Int,Int}(6=>1)

    context("updating best candidate index if old and new non-dominated") do
        a = EpsBoxArchive(scheme, max_size=100)

        BlackBoxOptim.add_candidate!(a, convert(IndexedTupleFitness, (1.0, 0.0), scheme), [0.0, 1.0], 1)
        @fact a.best_candidate_ix --> 1
        @fact length(a)           --> 1
        BlackBoxOptim.add_candidate!(a, convert(IndexedTupleFitness, (0.6, 0.6), scheme), [0.0, 2.0], 2)
        @fact a.best_candidate_ix --> 1
        @fact length(a)           --> 2
        BlackBoxOptim.add_candidate!(a, convert(IndexedTupleFitness, (0.0, 1.1), scheme), [0.0, 3.0], 3)
        @fact a.best_candidate_ix --> 1
        @fact length(a)           --> 3
        BlackBoxOptim.add_candidate!(a, convert(IndexedTupleFitness, (0.0, 0.9), scheme), [0.0, 4.0], 4)
        @fact a.best_candidate_ix --> 3
        @fact length(a)           --> 3
    end
  end

  context("shaffer1() Pareto frontier") do
      schaffer1(x) = (sumabs2(x), sumabs2(x .- 2.0))
      scheme = EpsBoxDominanceFitnessScheme{2}(0.1, is_minimizing=true)
      a = EpsBoxArchive(scheme)

      for t in linspace(0, 2, 1000)
          params = [t, t]
          BlackBoxOptim.add_candidate!(a, schaffer1(params), params)
      end

      @fact length(a) --> 40
  end

  context("archive copies the individuals") do
    a = EpsBoxArchive(scheme, max_size=100)
    @fact length(a) --> 0

    indiv = [1.5, 2.0]
    BlackBoxOptim.add_candidate!(a, (1.5, 2.0), indiv)
    # equal but not identical
    @fact a.frontier[1].params == indiv --> true
    @fact a.frontier[1].params === indiv --> false

    # modify the vector
    indiv[1] = 5.0
    # test that archived version is still the same
    @fact a.frontier[1].params --> [1.5, 2.0]
  end

end
