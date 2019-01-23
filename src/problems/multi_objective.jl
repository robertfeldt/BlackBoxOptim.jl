"""
``(N-1)``-dimensional manifold in ``N``-dimensional space.
"""
struct Hypersurface{N,SS<:SearchSpace}
    manifold::Function
    parameter_space::SS

    function Hypersurface(manifold::Function, parameter_space::SS) where {SS<:SearchSpace}
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
@generated function generate(surf::Hypersurface{N},
                             fs::EpsBoxDominanceFitnessScheme{N,F},
                             ::Type{Val{NP}},
                             param_step::Vector{F} = 0.1*fs.ϵ) where {N,F,NP}
    quote
        #archive = EpsBoxArchive(fs)
        pf = Dict{NTuple{N,Int}, IndexedTupleFitness{N,F}}()
        param = fill!(Vector{Float64}(undef, N-1), 0.0)
        #hat_compare = HatCompare(fs)
        Base.Cartesian.@nloops $(N-1) t d->range(surf.parameter_space.mins[d], stop=surf.parameter_space.maxs[d],
                                                 length=ceil(Int, surf.parameter_space.deltas[d]/param_step[d])) d->param[d]=t_d begin
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
function nondominated(fitnesses, fit_scheme::EpsBoxDominanceFitnessScheme{N,F}) where {N,F}
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
@generated function distance(a::NTuple{N,F}, b::NTuple{N,F}) where {N,F}
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
IGD(A::Vector{NTuple{N,F}}, B::Vector{NTuple{N,F}}) where {N,F<:Number} =
    sum(a -> minimum(b -> distance(a, b), B), A) / length(A)

"""
    IGD(ref::Hypersurface, sol::Vector{FitIndividual}, [two_sided=true])

Average Euclidean distance from the exact Pareto frontier of the problem (`ref`)
to the solution (`sol`) produced by the optimization method.
If `two_sided` is on, returns the maximum of `IGD(sol, ref)` and `IGD(nondominated(ref), sol)`.
"""
function IGD(ref::Hypersurface{N}, sol::AbstractVector{T},
             fit_scheme::EpsBoxDominanceFitnessScheme{N,F},
             ::Type{Val{NP}},
             param_step::Vector{F} = 0.1*fit_scheme.ϵ, two_sided::Bool=true
) where {N,F<:Number,T<:FrontierIndividualWrapper,NP}
    @assert fitness_type(sol[1]) == NTuple{N,F}
    ref_indexed_pts = values(generate(ref, fit_scheme, Val{NP}, param_step))
    nondom_ref_pts = NTuple{N,F}[val.orig for val in nondominated(ref_indexed_pts, fit_scheme)]
    orig_sol_pts = NTuple{N,F}[fitness(candi) for candi in sol]
    # how much `sol` should be inflated to contain all nondom_ref_pts
    res = IGD(nondom_ref_pts, orig_sol_pts)
    if two_sided
        # how much `ref` should be inflated to contain all `sol` points
        return max(res, IGD(orig_sol_pts, NTuple{N,F}[val.orig for val in ref_indexed_pts]))
    else
        return res
    end
end

"""
CEC09 Unconstrained Problem 8 objective function.
"""
function CEC09_UP8(x::Vector{Float64})
    N = length(x)
    N >= 3 || throw(ArgumentError("CEC09_UP8(): insufficient dimensions"))
    sumix = 3
    nterms = zeros(Int, 3)
    sums = zeros(Float64, 3)
    for j in 3:N
        sums[sumix] += (x[j] - 2*x[2]*sin(2π*x[1] + π/N*j))^2
        nterms[sumix] += 1
        sumix += 1
        if sumix > 3
            sumix = 1
        end
    end
    for i in 1:3
        if nterms[i] > 0
            sums[i] *= 2/nterms[i]
        end
    end
    (cos(0.5*π*x[1]) * cos(0.5*π*x[2]) + sums[1],
     cos(0.5*π*x[1]) * sin(0.5*π*x[2]) + sums[2],
     sin(0.5*π*x[1]) + sums[3])
end

"""
CEC09 Unconstrained Problem 8 Pareto Frontier.
Parameterized by t[1]=f[1] and t[2]=f[2].
"""
function CEC09_UP8_PF(t::Vector{Float64}, ::Type{Val{NP}}) where NP
    d=sum(abs2, t)
    (t[1], t[2], d <= 1.0 ? (1.0 - d)^(1/3) : NaN)
end

"""
The collection of CEC09 unconstrained
multi-objective problems.

See http://dces.essex.ac.uk/staff/zhang/MOEAcompetition/cec09testproblem0904.pdf.pdf
"""
const CEC09_Unconstrained_Set = Dict{Int,FunctionBasedProblemFamily}(
    8 => FunctionBasedProblemFamily(CEC09_UP8, "CEC09 UP8",  ParetoFitnessScheme{3}(is_minimizing=true),
                                    (-2.0, 2.0),
                                    Hypersurface(CEC09_UP8_PF, symmetric_search_space(2, (0.0, 1.0))),
                                    symmetric_search_space(2, (0.0, 1.0)))
)

schaffer1(x) = (sum(abs2, x), sum(xx -> abs2(xx - 2.0), x))
schaffer1_PF(t, ::Type{Val{NP}}) where {NP} = (NP*t[1]^2, NP*(2-t[1])^2)

const Schaffer1Family = FunctionBasedProblemFamily(schaffer1, "Schaffer1", ParetoFitnessScheme{2}(is_minimizing=true),
                                (-10.0, 10.0),
                                Hypersurface(schaffer1_PF, symmetric_search_space(1, (0.0, 2.0))))
