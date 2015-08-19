include("bimodal_cauchy_distribution.jl")

const ADE_DefaultOptions = chain(DE_DefaultOptions, @compat Dict{Symbol,Any}(
  # Distributions we will use to generate new F and CR values.
  :fdistr => bimodal_cauchy(0.65, 0.1, 1.0, 0.1),
  :crdistr => bimodal_cauchy(0.1, 0.1, 0.95, 0.1),
  :SearchSpace => symmetric_search_space(1)
))

# Specific data and functions for adaptation
# An Adaptive DE typically changes parameters of the search dynamically. This is
# typically done in the tell! function when we know if the trial vector
# was better than the target vector.
type AdaptiveDiffEvoParameters
  # Distributions we will use to generate new F and CR values.
  # FIXME allow any distribution?
  fdistr::BimodalCauchy
  crdistr::BimodalCauchy

  fs::Vector{Float64}   # One f value per individual in population
  crs::Vector{Float64}  # One cr value per individual in population

  function AdaptiveDiffEvoParameters(fdistr::BimodalCauchy = bimodal_cauchy(0.65, 0.1, 1.0, 0.1),
                                     crdistr::BimodalCauchy = bimodal_cauchy(0.1, 0.1, 0.95, 0.1))
    new(fdistr, crdistr,
        Array(Float64,0), Array(Float64,0)) # start with empty arrays because the population size unknown
  end
end

AdaptiveDiffEvoParameters(options::Parameters) = AdaptiveDiffEvoParameters(options[:fdistr], options[:crdistr])

function crossover_parameters(params::AdaptiveDiffEvoParameters, pop::Population, target_index)
  # initialize the f & cr array
  if isempty(params.fs)
    params.fs = [sample_bimodal_cauchy(params.fdistr; truncateBelow0 = false) for i in 1:popsize(pop)]
  end
  if isempty(params.crs)
    params.crs = [sample_bimodal_cauchy(params.crdistr; truncateBelow0 = false) for i in 1:popsize(pop)]
  end
  return (params.crs[target_index], params.fs[target_index])
end

function adjust!(params::AdaptiveDiffEvoParameters, index, is_improved::Bool)
    if !is_improved
      # The trial vector for this target was not better so we change the f and cr constants.
      params.fs[index] = sample_bimodal_cauchy(params.fdistr; truncateBelow0 = false)
      params.crs[index] = sample_bimodal_cauchy(params.crdistr; truncateBelow0 = false)
    end
end

# An Adaptive DE crossover operator changes cr and f parameters of the search dynamically.
type AdaptiveDiffEvoRandBin{N} <: DiffEvoCrossoverOperator{N,1}
  params::AdaptiveDiffEvoParameters

  AdaptiveDiffEvoRandBin(params::AdaptiveDiffEvoParameters) = new(params)
end

Base.convert{N}(::Type{AdaptiveDiffEvoRandBin{N}}, options::Parameters) =
    AdaptiveDiffEvoRandBin{N}(AdaptiveDiffEvoParameters(options))

crossover_parameters(xover::AdaptiveDiffEvoRandBin, pop::Population, target_index) =
    crossover_parameters(xover.params, pop, target_index)

adjust!{F}(xover::AdaptiveDiffEvoRandBin, op_index::Int, candi_index::Int,
                 new_fitness::F, old_fitness::F, is_improved::Bool) =
    adjust!(xover.params, candi_index, is_improved)

typealias AdaptiveDiffEvoRandBin1 AdaptiveDiffEvoRandBin{3}
typealias AdaptiveDiffEvoRandBin2 AdaptiveDiffEvoRandBin{5}

function adaptive_diffevo(problem::OptimizationProblem,
                 options::Parameters, name::ASCIIString,
                 select::IndividualsSelector = SimpleSelector(),
                 crossover::DiffEvoCrossoverOperator =
                    convert(AdaptiveDiffEvoRandBin1, chain(ADE_DefaultOptions, options)))
  opts = chain(ADE_DefaultOptions, options)
  pop = population(problem, opts)
  DiffEvoOpt(name, pop, select, crossover,
             RandomBound(search_space(problem)))
end

adaptive_de_rand_1_bin(problem::OptimizationProblem, options::Parameters = EMPTY_PARAMS,
              name = "AdaptiveDE/rand/1/bin") =
    adaptive_diffevo(problem, options, name)

adaptive_de_rand_1_bin_radiuslimited(problem::OptimizationProblem, options::Parameters = EMPTY_PARAMS,
                                     name = "AdaptiveDE/rand/1/bin/radiuslimited") =
    adaptive_diffevo(problem, options, name,
                     RadiusLimitedSelector(chain(ADE_DefaultOptions, options)[:SamplerRadius]))
