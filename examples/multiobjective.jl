using BlackBoxOptim, Gadfly

# run Borg MOAE
schaffer1(x) = (sumabs2(x), sumabs2(x .- 2.0))
res = bboptimize(schaffer1; Method=:borg_moea,
                FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true),
                SearchRange=(-10.0, 10.0), NumDimensions=3, ϵ=0.1,
                MaxSteps=15000, TraceInterval=1.0, TraceMode=:verbose);

# the parameterized version of the exact Pareto frontier
# N is the number of the problem (not fitness) dimensions
pareto_curve_func(t, ::Type{Val{N}}) where {N} = (N*t[1]^2, N*(2-t[1])^2)
pareto_curve = BlackBoxOptim.Hypersurface(pareto_curve_func,
                                          RectSearchSpace(1, (0.0, 2.0)))

# generate the set of ϵ-indexed points on the exact Pareto frontier
pareto_pts = BlackBoxOptim.generate(pareto_curve,
                                   fitness_scheme(res), Val{numdims(res)});
# calculate the distance between the solution and the exact frontier
BlackBoxOptim.IGD(pareto_curve, pareto_frontier(res), fitness_scheme(res), Val{numdims(res)})

# draw the results
plot(layer(x = [x.orig[1] for x in values(pareto_pts)], # Draw the exact Pareto frontier
           y = [x.orig[2] for x in values(pareto_pts)],
           Geom.line),
    layer(x= [fitness(x).orig[1] for x in pareto_frontier(res)], # the result of the method
          y= [fitness(x).orig[2] for x in pareto_frontier(res)],
          color= [haskey(pareto_pts, fitness(x).index) for x in pareto_frontier(res)],
          Geom.point))
