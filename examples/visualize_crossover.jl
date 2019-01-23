using BlackBoxOptim, Plots, Random

"""
    plot_crossover(xover::CrossoverOperator, parents::Matrix{Float64};
                   n=5000, shuffle_parents=false)

Plot the distribution of offsprings for the crossover operator `xover`
and given `parents`.
"""
function plot_crossover(xover::CrossoverOperator, parents::Matrix{Float64};
                        n=5000, shuffle_parents=false)
    @assert size(parents, 2) == numparents(xover)
    parents_pop = FitPopulation(parents, NaN)
    parent_ixs = collect(1:numparents(xover))
    children_mtx = hcat([hcat(apply!(xover, [fill!(Individual(undef, size(parents, 1)), NaN) for _ in 1:numchildren(xover)],
                              zeros(Int, numchildren(xover)), parents_pop,
                              shuffle_parents ? shuffle(parent_ixs) : parent_ixs)...) for _ in 1:n]...)
    plot(view(children_mtx, 1, :), view(children_mtx, 2, :),
         seriestype=:scatter, title=string(typeof(xover)), label="children",
         markersize=1.5, markerstrokewidth=0, markercolor=colorant"orange")
    scatter!(view(parents_pop.individuals, 1, :), view(parents_pop.individuals, 2, :),
             markersize=5, markercolor=colorant"purple", label="parents")
end

plot_crossover(SimplexCrossover{3}(1.2),
               [[0.5, 0.0] [-0.2, 0.0] [0.0, 1.0]])

plot_crossover(SimulatedBinaryCrossover(0.2, 5.0),
               [[1.0, 0.2] [-0.2, 1.0]])

plot_crossover(UnimodalNormalDistributionCrossover{2}(0.5, 0.05),
               [[0.5, 0.0] [1.0, 3.0]])

plot_crossover(UnimodalNormalDistributionCrossover{3}(0.5, 0.05),
               [[0.5, 0.0] [1.0, 3.0] [3.0, 1.0]])

plot_crossover(ParentCentricCrossover{3}(0.1, 0.2),
               [[0.5, 0.0] [1.0, 3.0] [3.0, 1.0]])

plot_crossover(ParentCentricCrossover{3}(0.1, 0.2),
               [[0.5, 0.0] [1.0, 3.0] [3.0, 1.0]], shuffle_parents=true)
