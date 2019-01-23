using BlackBoxOptim: fitness, index, tag, rank_by_fitness!, bbsetup
using BlackBoxOptim: hassolution, setsolution!, solution

# The ask tell interface aims to allow flexibility in cases where
# the evaluator might not be easy/possible to define up-front
# or when the optimization needs to be "driven" from the outside
# rather than from BBO itself.
@testset "Top-level ask/tell interface" begin
    @testset "Single float fitness, with intermediate solutions" begin

        # A typical use case for this is that we need to map parameter sets to
        # intermediate "solutions" before we can evaluate their fitness. We want
        # to avoid the repeated generation of solutions since that might be costly.
        # In this example it is not but in general it might be, and then we want
        # to "cache" the intermediate solutions with the parameters.
        params2solution(params) = string(sum(params))
        fitn(sol) = length(sol)

        # Since the user will step through the ask/tell cycle themselves we should
        # not require knowledge of the optimization function / evaluator. For cases
        # where the fitness is not a single Float64 one should state the fitness
        # scheme explicitly as for the multi-objective case?
        oc = bbsetup(; SearchRange = (-1.0, 1.0), NumDimensions = 5)

        # When we ask for solutions we get either solutions for params
        # we previously supplied or empty Nullables. The alternative would
        # be to embed the solutions in the Candidate struct and let user
        # have getters/setters for them. Seems unnecessary to add this 
        # level of complexity. Easier to have a general top-level interface
        # where we can get and set solutions/phenotypes associated with each 
        # set of parameters. OTOH we will probably save the solutions in the
        # candidates anyway so users that require a more control can just ask
        # to get the candidates themselves and then return them (sorted).
        # Here we show the simpler, first version. It will unpack the params
        # and the solutions from the (internal) candidates and thus does not
        # require the user to know anything about the Candidate type or internals
        # of BBO.
        params, solutions = ask(oc; withSolutions = true)

        @test length(params) > 0
        @test length(params) == length(solutions)

        # In this case there are no solutions yet but in general there might be.
        # Our task is now to sort the params and solutions based on their fitness.
        fs = Array{Float64}()
        for i in eachindex(params)
            p = params[i]
            s = solutions[i]
            @test typeof(p) <: Vector{Float64}
            @test typeof(s) <: Nullable
            @test fitness(c) == NaN       # Since not yet evaluated (and since we did not specify a non-standard fitness scheme).
            @test isnull(s)               # Since not yet evaluated. If one of the params is one for which we have previously supplied the solution is should be returned (in a Nullable) instead.

            ns = isnull(s) ? params2solution(p) : s # Don't calc new solution is already given before.
            push!(newsolutions, ns)
            push!(fs, fitn(ns))
        end

        # A usability risk here is that the user only sorts one of the arrays
        # returned? Another risk is that the fitnesses might not be comparable on 
        # a global scale, i.e. between calls to ask/tell. This is likely to mess up
        # tracing and archives as well as the saving of the "best" individuals. OTOH,
        # since optimization is "driven" from the outside maybe the resposibility for 
        # keeping history, tracing etc should be on the user rather than on BBO?
        p = sortperm(fs)
        tell!(oc, params[p], fs[p], newsolutions[p])
    end
end