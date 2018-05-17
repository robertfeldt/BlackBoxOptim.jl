"""
Single objective optimization methods accepted by `bboptimize()`.

The values are the method initialization routines or types derived from `Optimizer`.
"""
const SingleObjectiveMethods = ParamsDict(
    :random_search => random_search,
    :de_rand_1_bin => de_rand_1_bin,
    :de_rand_2_bin => de_rand_2_bin,
    :de_rand_1_bin_radiuslimited => de_rand_1_bin_radiuslimited,
    :de_rand_2_bin_radiuslimited => de_rand_2_bin_radiuslimited,
    :adaptive_de_rand_1_bin => adaptive_de_rand_1_bin,
    :adaptive_de_rand_1_bin_radiuslimited => adaptive_de_rand_1_bin_radiuslimited,
    :separable_nes => separable_nes,
    :xnes => xnes,
    :dxnes => dxnes,
    :resampling_memetic_search => resampling_memetic_searcher,
    :resampling_inheritance_memetic_search => resampling_inheritance_memetic_searcher,
    :simultaneous_perturbation_stochastic_approximation => SimultaneousPerturbationSA2,
    :generating_set_search => GeneratingSetSearcher,
    :probabilistic_descent => direct_search_probabilistic_descent,
)

const SingleObjectiveMethodNames = sort!(collect(keys(SingleObjectiveMethods)))

"""
Multi-objective optimization methods accepted by `bboptimize()`.

The values are the method initialization routines or types derived from `Optimizer`.
"""
const MultiObjectiveMethods = ParamsDict(
    :borg_moea => borg_moea
)

const MultiObjectiveMethodNames = sort!(collect(keys(MultiObjectiveMethods)))

"""
Names of optimization methods accepted by `bboptimize()`,
`:Method` keyword argument.
"""
const MethodNames = sort!(vcat(SingleObjectiveMethodNames,
                               MultiObjectiveMethodNames))
