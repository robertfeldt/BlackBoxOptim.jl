facts("EpsBoxArchive") do
    # check that frontier elements are mutually nondominated and
    # ther are no epsilon-index duplicates
    function check_frontier(a::EpsBoxArchive)
        for (i, el_i) in enumerate(pareto_frontier(a))
            for (j, el_j) in enumerate(pareto_frontier(a))
                if j > i
                    comp, index_match = hat_compare(fitness(el_i), fitness(el_j), fitness_scheme(a))
                    if comp != 0 || index_match
                        return comp, index_match # return if violating pair is found
                    end
                end
            end
        end
        return (0, false) # all non-dominated, no index matches
    end

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
    @fact BlackBoxOptim.noprogress_streak(a) --> 0
    @fact BlackBoxOptim.tagcounts(a) --> Dict{Int,Int}()

    BlackBoxOptim.add_candidate!(a, convert(IndexedTupleFitness, (2.21, 1.12), scheme), [0.0, 5.0], 1, 8)
    @fact capacity(a)         --> 100
    @fact length(a)           --> 1
    @fact a.best_candidate_ix --> 1
    @fact a.frontier[1].num_fevals --> 8
    @fact BlackBoxOptim.noprogress_streak(a) --> 0
    @fact BlackBoxOptim.tagcounts(a) --> Dict{Int,Int}(1=>1)

    BlackBoxOptim.add_candidate!(a, (3.21, 1.12), [1.0, 5.0], 2)
    @fact capacity(a)         --> 100
    @fact length(a)           --> 1
    @fact a.best_candidate_ix --> 1
    @fact a.frontier[1].num_fevals --> 8
    @fact BlackBoxOptim.noprogress_streak(a) --> 1
    @fact BlackBoxOptim.tagcounts(a) --> Dict{Int,Int}(1=>1)

    BlackBoxOptim.add_candidate!(a, (1.21, 2.11), [1.2, 3.0], 3)
    @fact capacity(a)         --> 100
    @fact length(a)           --> 2
    @fact a.frontier[2].num_fevals --> 3
    @fact a.best_candidate_ix --> 2
    @fact BlackBoxOptim.noprogress_streak(a) --> 0
    @fact BlackBoxOptim.tagcounts(a) --> Dict{Int,Int}(1=>1, 3=>1)

    BlackBoxOptim.add_candidate!(a, (0.52, 3.15), [1.5, 3.0], 3)
    @fact capacity(a)         --> 100
    @fact length(a)           --> 3
    @fact a.best_candidate_ix --> 2
    @fact BlackBoxOptim.noprogress_streak(a) --> 0
    @fact BlackBoxOptim.tagcounts(a) --> Dict{Int,Int}(1=>1, 3=>2)

    BlackBoxOptim.add_candidate!(a, (0.22, 0.22), [1.5, 3.0], 5)
    @fact capacity(a)         --> 100
    @fact length(a)           --> 1
    @fact a.best_candidate_ix --> 1
    @fact BlackBoxOptim.noprogress_streak(a) --> 0
    @fact BlackBoxOptim.tagcounts(a) --> Dict{Int,Int}(5=>1)

    BlackBoxOptim.add_candidate!(a, (0.21, 0.21), [1.5, 3.0], 6)
    @fact capacity(a)         --> 100
    @fact length(a)           --> 1
    @fact a.best_candidate_ix --> 1
    @fact BlackBoxOptim.noprogress_streak(a) --> 1
    @fact BlackBoxOptim.tagcounts(a) --> Dict{Int,Int}(6=>1)

    context("noprogress_streak(since_restart=0) is reset after restart") do
        @fact BlackBoxOptim.noprogress_streak(a, since_restart=false) --> 1
        @fact BlackBoxOptim.noprogress_streak(a, since_restart=true) --> 1
        BlackBoxOptim.notify!(a, :restart)
        @fact BlackBoxOptim.noprogress_streak(a, since_restart=false) --> 1
        @fact BlackBoxOptim.noprogress_streak(a, since_restart=true) --> 0
        BlackBoxOptim.add_candidate!(a, (0.2, 0.2), [1.5, 3.0], 6)
        @fact BlackBoxOptim.noprogress_streak(a, since_restart=false) --> 2
        @fact BlackBoxOptim.noprogress_streak(a, since_restart=true) --> 1
    end

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

    context("handling dulicate elements") do
        a = EpsBoxArchive(scheme, max_size=100)

        BlackBoxOptim.add_candidate!(a, convert(IndexedTupleFitness, (1.25, 0.0), scheme), [0.0, 1.0], 1)
        @fact a.best_candidate_ix --> 1
        @fact length(a)           --> 1
        BlackBoxOptim.add_candidate!(a, convert(IndexedTupleFitness, (0.6, 0.6), scheme), [0.0, 2.0], 2)
        @fact a.best_candidate_ix --> 2
        @fact length(a)           --> 2
        BlackBoxOptim.add_candidate!(a, convert(IndexedTupleFitness, (0.0, 1.3), scheme), [0.0, 3.0], 3)
        @fact a.best_candidate_ix --> 2
        @fact length(a)           --> 3
        BlackBoxOptim.add_candidate!(a, convert(IndexedTupleFitness, (0.6, 0.6), scheme), [0.0, 4.0], 4)
        @fact a.best_candidate_ix --> 2
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

  context("simulate optimization") do
    # generate random fitness from (0..2, 0..2, 0..2)
    # and lying inside the (1.0-r1)^1.5 + (1.0-r2)^1.5 + (1.0-r3)^1.5 <= 8.0
    function rfitness()
        while true
            r1 = rand()
            r2 = rand()
            r3 = rand()
            if (1.0-r1)^1.5 + (1.0-r2)^1.5 + (1.0-r3)^1.5 <= 1.0
                return (2.0*r1, 2.0*r2, 2.0*r3)
            end
        end
    end

    scheme3 = EpsBoxDominanceFitnessScheme{3}(0.1, is_minimizing=true)
    a = EpsBoxArchive(scheme3, max_size=1000)
    indiv = [0.0, 0.0]
    for i in 1:10000
        BlackBoxOptim.add_candidate!(a, rfitness(), indiv)
        # periodically check that pareto_frontier(a) is indeed
        # epsilon-index frontier
        if mod(i, 100) == 0
            @fact length(a) --> isposdef
            @fact check_frontier(a) --> (0, false)
        end
    end
  end
end
