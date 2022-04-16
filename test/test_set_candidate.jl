using BlackBoxOptim: set_candidate!, candidate, set_candidates!

# a 2d optimization problem with an optimum at (3.14, 7.2)
fixed_optimum_prob(x) = (x[1] - 3.14)^2 + (x[2] - 7.2)^4
const FixedOptimum = [3.14, 7.2]
const FitnessOptimum = fixed_optimum_prob(FixedOptimum)

@testset "set_candidate! for DE population optimizers" begin
    for m in [:de_rand_1_bin, :de_rand_2_bin, :de_rand_1_bin_radiuslimited, :de_rand_2_bin_radiuslimited, 
                :adaptive_de_rand_1_bin, :adaptive_de_rand_1_bin_radiuslimited]
        PopSize = 10
        b = bbsetup(fixed_optimum_prob; Method = m, MaxFuncEvals = 10*PopSize, # Give it a chance to be sampled so best fitness is in archive
            NumDimensions = 2, SearchRange = (-10.0, 10.0), PopulationSize = PopSize)
        set_candidate!(b.optimizer, FixedOptimum)
        inds = population(b.optimizer).individuals
        
        # There is at least one position in the population that is set to the FixedOptimum:
        @test any(i -> inds[:, i] == FixedOptimum, 1:PopSize)

        res = bboptimize(b; TraceMode = :silent)
        @test isapprox(best_fitness(res), FitnessOptimum)
    end
end

@testset "set_candidate! for stepping optimizers" begin
    for m in [:generating_set_search, :resampling_memetic_search, :resampling_inheritance_memetic_search,
                :simultaneous_perturbation_stochastic_approximation]
        b = bbsetup(fixed_optimum_prob; Method = m, MaxFuncEvals = 2, 
            NumDimensions = 2, SearchRange = (-10.0, 10.0), PopulationSize = 5)
        set_candidate!(b.optimizer, FixedOptimum)
        @test candidate(b.optimizer) == FixedOptimum

        # Note that we exclude simultaneous_perturbation_stochastic_approximation since it updates even
        # if new candidates give worse fitness => not guaranteed to find the optimum. We should debug this at some point.
        if m != :simultaneous_perturbation_stochastic_approximation
            res = bboptimize(b; TraceMode = :silent)
            @test isapprox(best_fitness(res), FitnessOptimum)
        end
    end
end

@testset "set_candidate! for NES optimizers" begin
    for m in [:dxnes, :xnes, :separable_nes]
        b = bbsetup(fixed_optimum_prob; Method = m, MaxFuncEvals = 100, 
            NumDimensions = 2, SearchRange = (-10.0, 10.0), PopulationSize = 5)
        set_candidate!(b.optimizer, FixedOptimum)
        @test candidate(b.optimizer) == FixedOptimum

        res = bboptimize(b; TraceMode = :silent)
        @test isapprox(best_fitness(res), FitnessOptimum, atol = 1e-1)
    end
end

@testset "set_candidate! for BORG" begin
    # Best aggregated fitness should be for [0.5, 0.5, 0.5]:
    fitness_2obj(x) = (sum(abs2, x), sum(abs2, x .- 1.0))
    x0 = 0.5 * ones(3)

    PopSize = 10
    b = bbsetup(fitness_2obj; Method=:borg_moea,
            FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true),
            MaxFuncEvals = 10*PopSize, PopulationSize = PopSize,
            SearchRange=(-10.0, 10.0), NumDimensions=3, ϵ=0.05);

    set_candidate!(b.optimizer, x0)

    inds = b.optimizer.population.individuals
    # There is at least one position in the population that is set to x0
    @test any(i -> inds[:, i] == x0, 1:PopSize)

    # Now ensure we actually get back (sum(abs2, x0), sum(abs2, x0 .- 1.0)) as the best fitness (since it is the best aggregated one).
    res = bboptimize(b; TraceMode = :silent)
    @test best_fitness(res) == fitness_2obj(x0) 
end

@testset "set_candidates! for DE population optimizers" begin
    for m in [:de_rand_1_bin, :de_rand_2_bin, :de_rand_1_bin_radiuslimited, :de_rand_2_bin_radiuslimited, 
                :adaptive_de_rand_1_bin, :adaptive_de_rand_1_bin_radiuslimited]
        PopSize = 10
        b = bbsetup(fixed_optimum_prob; Method = m, MaxFuncEvals = 10*PopSize, # Give it a chance to be sampled so best fitness is in archive
            NumDimensions = 2, SearchRange = (-10.0, 10.0), PopulationSize = PopSize)

        # Get a random population and then add the optimum to ensure it is in there
        randval() = -10.0 + rand() * 20.0
        initial_population = [[randval(), randval()] for _ in 1:(PopSize-rand(1:2))]
        push!(initial_population, FixedOptimum)

        # Now set the population
        set_candidates!(b.optimizer, initial_population)

        # There is at least one position in the population that is set to the FixedOptimum:
        inds = population(b.optimizer).individuals
        @test any(i -> inds[:, i] == FixedOptimum, 1:PopSize)

        # and all of the starting points can be found somewhere in the population
        for startingpoint in initial_population
            @test any(i -> inds[:, i] == startingpoint, 1:PopSize)
        end
    
        res = bboptimize(b; TraceMode = :silent)
        @test isapprox(best_fitness(res), FitnessOptimum)
    end
end

@testset "set_candidates! when calling bboptimize directly" begin
    PopSize = 10

    for m in [:de_rand_1_bin, :de_rand_2_bin, :de_rand_1_bin_radiuslimited, :de_rand_2_bin_radiuslimited, 
        :adaptive_de_rand_1_bin, :adaptive_de_rand_1_bin_radiuslimited]
        
        # Get a random population and then add the optimum to ensure it is in there
        randval() = -10.0 + rand() * 20.0
        startingpoints = [[randval(), randval()] for _ in 1:(PopSize-rand(1:2))]
        push!(startingpoints, FixedOptimum)

        res = bboptimize(fixed_optimum_prob, startingpoints; Method = m, MaxFuncEvals = 10*PopSize, # Give it a chance to be sampled so best fitness is in archive
            NumDimensions = 2, SearchRange = (-10.0, 10.0), PopulationSize = PopSize,
            TraceMode = :silent)
        @test isapprox(best_fitness(res), FitnessOptimum)
    end
end