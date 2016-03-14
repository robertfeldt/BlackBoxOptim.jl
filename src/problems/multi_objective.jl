"""
  ``(N-1)``-dimensional manifold in ``N``-dimensional space.
"""
immutable Hypersurface{N,SS<:SearchSpace}
    manifold::Function
    parameter_space::SS

    function Base.call{SS<:SearchSpace}(::Type{Hypersurface}, manifold::Function, parameter_space::SS)
        N = numdims(parameter_space)+1
        int_pt = manifold(0.5*(mins(parameter_space)+maxs(parameter_space)), Val{1})
        length(int_pt) == N || throw(DimensionMismatch())
        new{N,SS}(manifold, parameter_space)
    end
end

"""
    generate(surf::Hypersurface, fs::EpsBoxDominanceFitnessScheme, param_step = 0.1*fs.ϵ)

    Generate the points of the hypersurface using the
    discretization defined by ϵ-box fitness schema.
"""
@generated function generate{N,F,NP}(surf::Hypersurface{N}, fs::EpsBoxDominanceFitnessScheme{N,F}, ::Type{Val{NP}}, param_step::Number = 0.1*fs.ϵ)
    quote
        #archive = EpsBoxArchive(fs)
        pf = Dict{NTuple{N,Int}, IndexedTupleFitness{N,F}}()
        param = Vector{Float64}(N-1)
        #hat_compare = HatCompare(fs)
        Base.Cartesian.@nloops $(N-1) t d->linspace(surf.parameter_space.mins[d], surf.parameter_space.maxs[d], ceil(Int, surf.parameter_space.deltas[d]/param_step)) d->param[d]=t_d begin
            fit = surf.manifold(param, Val{NP})
            if !isnafitness(fit, fs) # NA if given parameters do not correspond to any point on the manifold
                ifit = convert(IndexedTupleFitness, fit, fs)
                pf[ifit.index] = ifit
                #add_candidate!(archive, convert(IndexedTupleFitness, fit, fs), Individual())
            end
        end
        return pf
    end
end

"""
    nondominated(fitnesses, fit_scheme)

    Filter `fitnesses` removing all dominated values.
"""
function nondominated{N,F}(fitnesses::Vector{IndexedTupleFitness{N,F}}, fit_scheme::EpsBoxDominanceFitnessScheme{N,F})
    arch = EpsBoxArchive(fit_scheme)
    res = sizehint!(Vector{IndexedTupleFitness{N,F}}(), length(fitnesses))
    empty_params = Individual()
    for f in fitnesses
        add_candidate!(arch, f, empty_params)
    end
    IndexedTupleFitness{N,F}[fitness(pt) for pt in pareto_frontier(arch)]
end

"""
    distance(a::NTuple{N,F}, b::NTuple{N,F})

    Euclidean distance from `a` to `b`.
"""
@generated function distance{N,F}(a::NTuple{N,F}, b::NTuple{N,F})
    quote
        res = 0.0
        Base.Cartesian.@nexprs $N i -> res += (a[i]-b[i])^2
        return sqrt(res)
    end
end

"""
    IGD(A::Vector{NTuple{N,F}}, B::Vector{NTuple{N,F}}, [two_sided=true])

    The average Euclidean distance from the points of `A` to the points of `B`.
"""
IGD{N,F<:Number}(A::Vector{NTuple{N,F}}, B::Vector{NTuple{N,F}}) =
    sum(a -> minimum(b -> distance(a, b), B), A) / length(A)

"""
    IGD(ref::Hypersurface, sol::Vector{EpsBoxFrontierIndividual}, [two_sided=true])

    Average Euclidean distance from the exact Pareto frontier of the problem (`ref`)
    to the solution (`sol`) produced by the optimization method.
    If `two_sided` is on, returns the maximum of `IGD(sol, ref)` and `IGD(nondominated(ref), sol)`.
"""
function IGD{N,F<:Number,NP}(ref::Hypersurface{N}, sol::Vector{EpsBoxFrontierIndividual{N,F}},
                             fit_scheme::EpsBoxDominanceFitnessScheme{N,F},
                             ::Type{Val{NP}},
                             param_step::Number = 0.1*fit_scheme.ϵ, two_sided::Bool=true)
    ref_pts = values(generate(ref, fit_scheme, Val{NP}, param_step))
    orig_sol_pts = NTuple{N,F}[fitness(candi).orig for candi in sol]
    nondom_ref_pts = nondominated(IndexedTupleFitness{N,F}[val for val in ref_pts], fit_scheme)
    # how much `sol` should be inflated to contain all nondom_ref_pts
    res = IGD(NTuple{N,F}[val.orig for val in nondom_ref_pts], orig_sol_pts)
    if two_sided
        # how much `ref` should be inflated to contain all `sol` points
        return max(res, IGD(orig_sol_pts, NTuple{N,F}[val.orig for val in ref_pts]))
    else
        return res
    end
end

schaffer1(x) = (sumabs2(x), sumabs2(x .- 2.0))
schaffer1_PF{NP}(t, ::Type{Val{NP}}) = (NP*t[1]^2, NP*(2-t[1])^2)

const Schaffer1Family = FunctionBasedProblemFamily(schaffer1, "Schaffer1", ParetoFitnessScheme{2}(is_minimizing=true),
                                (-10.0, 10.0),
                                Hypersurface(schaffer1_PF, symmetric_search_space(1, (0.0, 2.0))))
