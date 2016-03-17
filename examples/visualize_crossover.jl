using BlackBoxOptim, Gadfly

"""
    plot_crossover(xover::CrossoverOperator, parents::Matrix{Float64};
                   n=5000, shuffle_parents=false)

    Plot the distribution of offsprings for the crossover operator `xover`
    and given `parents`.
"""
function plot_crossover{NP,NC}(xover::CrossoverOperator{NP,NC}, parents::Matrix{Float64};
                        n=5000, shuffle_parents=false)
    @assert size(parents, 2) == NP
    parents_pop = FitPopulation(parents, NaN)
    parent_ixs = collect(1:NP)
    children_mtx = hcat([hcat(apply!(xover, [Individual(size(parents, 1)) for _ in 1:NC],
                              zeros(Int, NC), parents_pop,
                              shuffle_parents ? shuffle(parent_ixs) : parent_ixs)...) for _ in 1:n]...)
    plot(layer(x=parents_pop.individuals[1,:], y=parents_pop.individuals[2,:], Geom.point, # parents
               Theme(default_color=colorant"purple", default_point_size=8.0pt, highlight_width=0pt)),
         layer(x=children_mtx[1,:], y=children_mtx[2,:], Geom.point,
               Theme(default_color=colorant"orange", default_point_size=1.5pt, highlight_width=0pt))
    )
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
