include("bimodal_cauchy_distribution.jl")

ADE_DefaultOptions = mergeparam(DE_DefaultOptions, @compat Dict{Symbol,Any}(
  # Distributions we will use to generate new F and CR values.
  :fdistr => bimodal_cauchy(0.65, 0.1, 1.0, 0.1),
  :crdistr => bimodal_cauchy(0.1, 0.1, 0.95, 0.1),
  :SearchSpace => symmetric_search_space(1)
))

# Specific data and functions for adaptation
# An Adaptive DE typically changes parameters of the search dynamically. This is
# typically done in the tell! function when we know if the trial vector
# was better than the target vector.
type AdaptiveDiffEvoParameters <: DiffEvoParameters
  # Distributions we will use to generate new F and CR values.
  # FIXME allow any distribution?
  fdistr::BimodalCauchy
  crdistr::BimodalCauchy

  fs::Vector{Float64}   # One f value per individual in population
  crs::Vector{Float64}  # One cr value per individual in population

  function AdaptiveDiffEvoParameters(popsize::Int = 100,
                                     fdistr::BimodalCauchy = bimodal_cauchy(0.65, 0.1, 1.0, 0.1),
                                     crdistr::BimodalCauchy = bimodal_cauchy(0.1, 0.1, 0.95, 0.1) )
    new(fdistr, crdistr,
        [sample_bimodal_cauchy(fdistr; truncateBelow0 = false) for i in 1:popsize],
        [sample_bimodal_cauchy(crdistr) for i in 1:popsize])
  end
end

AdaptiveDiffEvoParameters(options, popsize::Int) = AdaptiveDiffEvoParameters(popsize, options[:fdistr], options[:crdistr])

crossover_parameters(params::AdaptiveDiffEvoParameters, index) = params.crs[index], params.fs[index]

function adjust!( params::AdaptiveDiffEvoParameters, index, is_improved::Bool )
    if !is_improved
      # The trial vector for this target was not better so we change the f and cr constants.
      params.fs[index] = sample_bimodal_cauchy( params.fdistr; truncateBelow0 = false)
      params.crs[index] = sample_bimodal_cauchy( params.crdistr; truncateBelow0 = false)
    end
end

function adaptive_diffevo(problem::OptimizationProblem, name::ASCIIString,
                 crossover::DiffEvoCrossoverOperator,
                 select::IndividualsSelector = SimpleSelector(),
                 options = @compat Dict{Symbol,Any}())
  opts = Parameters(options, ADE_DefaultOptions)
  pop = opts[:Population]
  DiffEvoOpt(name, pop, AdaptiveDiffEvoParameters(opts, popsize(pop)), select,
        NoMutation(), crossover, RandomBound(search_space(problem)))
end

adaptive_de_rand_1_bin(problem::OptimizationProblem, options = @compat(Dict{Symbol,Any}()),
              name = "AdaptiveDE/rand/1/bin") =
    adaptive_diffevo(problem, name, DiffEvoRandBin1(), SimpleSelector(), options)

adaptive_de_rand_1_bin_radiuslimited(problem::OptimizationProblem, options = @compat(Dict{Symbol,Any}()),
                                     name = "AdaptiveDE/rand/1/bin/radiuslimited") =
    adaptive_diffevo(problem, name, DiffEvoRandBin1(),
                     RadiusLimitedSelector(Parameters(options, ADE_DefaultOptions)[:SamplerRadius]), options)
