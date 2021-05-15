using BlackBoxOptim: set_candidate!, candidate

# a 2d optimization problem with an optimum at (3.14, 7.2)
fixed_optimum_prob(x) = (x[1] - 3.14)^2 + (x[2] - 7.2)^4
const FixedOptimum = [3.14, 7.2]
const FitnessOptimum = fixed_optimum_prob(FixedOptimum)

@testset "set_candidate! for DE population optimizers" begin
    for m in [:de_rand_1_bin, :de_rand_2_bin, :de_rand_1_bin_radiuslimited, :de_rand_2_bin_radiuslimited, 
                :adaptive_de_rand_1_bin, :adaptive_de_rand_1_bin_radiuslimited]
        println("Testing set_candidate! on $m")
        PopSize = 10
        b = bbsetup(fixed_optimum_prob; Method = m, MaxFuncEvals = 10*PopSize, # Give it a chance to be sampled so best fitness is in archive
            NumDimensions = 2, SearchRange = (-10.0, 10.0), PopulationSize = PopSize)
        set_candidate!(b.optimizer, FixedOptimum)
        inds = population(b.optimizer).individuals
        # There is at least one position in the population that is set to the FixedOptimum:
        @test any(i -> inds[:, i] == FixedOptimum, 1:PopSize)

        res = bboptimize(b)
        @test isapprox(best_fitness(res), FitnessOptimum)
    end
end

@testset "set_candidate! for stepping optimizers" begin
    for m in [:generating_set_search, :resampling_memetic_search, :resampling_inheritance_memetic_search,
                :simultaneous_perturbation_stochastic_approximation]
        println("Testing set_candidate! on $m")
        b = bbsetup(fixed_optimum_prob; Method = m, MaxFuncEvals = 2, 
            NumDimensions = 2, SearchRange = (-10.0, 10.0), PopulationSize = 5)
        set_candidate!(b.optimizer, FixedOptimum)
        @test candidate(b.optimizer) == FixedOptimum

        # Note that we exclude simultaneous_perturbation_stochastic_approximation since it updates even
        # if new candidates give worse fitness => not guaranteed to find the optimum. We should debug this at some point.
        if m != :simultaneous_perturbation_stochastic_approximation
            res = bboptimize(b)
            @test isapprox(best_fitness(res), FitnessOptimum)
        end
    end
end

@testset "set_candidate! for NES optimizers" begin
    for m in [:dxnes, :xnes, :separable_nes]
        println("Testing set_candidate! on $m")
        b = bbsetup(fixed_optimum_prob; Method = m, MaxFuncEvals = 100, 
            NumDimensions = 2, SearchRange = (-10.0, 10.0), PopulationSize = 5)
        set_candidate!(b.optimizer, FixedOptimum)
        @test candidate(b.optimizer) == FixedOptimum

        res = bboptimize(b)
        @test isapprox(best_fitness(res), FitnessOptimum, atol = 1e-1)
    end
end
