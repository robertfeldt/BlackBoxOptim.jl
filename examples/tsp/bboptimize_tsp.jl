using BlackBoxOptim

# Let's try to solve some TSP problems with continuous optimization
# algorithms by mapping from the float candidate vector to a TSP
# solution.
include("read_tsplib_file.jl")

TspProblemFile = "dantzig42.tsp" # Best known solution is 699
TspProblemFile = "fri26.tsp" # Best known solution is 937
const TSP = read_tsplib_file(TspProblemFile)

# We get inspiration from the mapping methods in section 2.4 of the paper:
# https://www.sciencedirect.com/science/article/pii/S2210650219304468?casa_token=Pyu6JJAzBQoAAAAA:J-EYa_kdwfDQidRAHbxmvnAwsmttUWO-hDESjTQyLlr70RMOfYdOijJs8Reao1yr-pXrIgpuxsY

abstract type FloatVectorMapping end
abstract type PermutationMapping <: FloatVectorMapping end

# Some mappings will use the best so far as a guide but by default not.
updatefitness!(p::PermutationMapping, fitness::Float64) = nothing

# This has also been called SPV (Smallest Position Value) but they
# look the same.
struct RankedOrderValueMapping <: PermutationMapping end

function apply(m::RankedOrderValueMapping, v::Vector{Float64})
    rov = zeros(Int, length(v))
    rov[sortperm(v)] = 1:length(v)
    rov
end
# v = [2.54, 1.52, 2.13, 3.48, 2.37, 2.87]
# @assert apply(RankedOrderValueMapping(), v) == [4, 1, 2, 6, 3, 5]

# sortperm right away should be the same though but faster:
struct SortPermMapping <: PermutationMapping
    reverse::Bool
    SortPermMapping(r::Bool = false) = new(r)
end
apply(m::SortPermMapping, v::Vector{Float64}) = sortperm(v, rev=m.reverse)

# The literature on mappings is strange with some papers
# basically just reversing the order and claiming novelty.
# For example: Li, Xiangtao, and Minghao Yin. "A hybrid cuckoo 
# search via Lévy flights for the permutation flow shop scheduling 
# problem." International Journal of Production Research 51.16 (2013): 4732-4754.
# Proposes Largest Ranked Value which just seems to be the reverse
# order of SPV/ROV above. I can't imagine this can have any effect
# on optimization performance.
struct LargestRankedValueMapping <: PermutationMapping end
function apply(m::LargestRankedValueMapping, v::Vector{Float64})
    rov = zeros(Int, length(v))
    rov[sortperm(v, rev=true)] = 1:length(v)
    rov
end

# Let's try something of our own, i.e. truncate/floor countinuous
# values to integers, then handle the duplicates by assigning the 
# remaining values based on which of the remaining ones it is closest
# to.
struct TruncateThenClosestRemainingMapping <: PermutationMapping
    N::Int
    s::Set{Int}
    TruncateThenClosestRemainingMapping(N::Int) = new(N, Set(1:N))
end
searchrange(m::TruncateThenClosestRemainingMapping) = (1.0, m.N + 1.0 - 1e-10)

function apply(m::TruncateThenClosestRemainingMapping, v::Vector{Float64})
    ints = floor.(Int, v)
    # Now repair by selecting the closest unselected
    selected = unique(ints)
    unselected = collect(setdiff(m.s, selected))
    counts = zeros(Int, m.N)
    for (idx, ival) in enumerate(ints)
        if counts[ival] > 0
            _, i = findmin(map(r -> abs(ival - r), unselected))
            r = unselected[i]
            counts[r] += 1
            ints[idx] = r
            deleteat!(unselected, i)
        else
            counts[ival] += 1
        end
    end
    ints
end

const ROV = RankedOrderValueMapping()
const SPM = SortPermMapping()
const LRV = LargestRankedValueMapping()
const TRC = TruncateThenClosestRemainingMapping(size(TSP))

function make_perm_mapping_fitness(tsp, pm::PermutationMapping)
    function tspfitnessfn(v::Vector{Float64})
        c = cost(tsp, apply(pm, v))
        updatefitness!(pm, c)
        c
    end
    tspfitnessfn
end

res_rov = bboptimize(make_perm_mapping_fitness(TSP, ROV); 
    SearchRange = (0.0, 10.0), NumDimensions = size(TSP),
    PopulationSize = 1000,
    MaxTime = 30.0)
# dantzig42.tsp best runs: 836, 818
# fri26.tsp best runs: 977

res_spm = bboptimize(make_perm_mapping_fitness(TSP, SPM);
    SearchRange = (0.0, 10.0), NumDimensions = size(TSP),
    PopulationSize = 1000,
    MaxTime = 30.0)
# dantzig42.tsp Best runs: 749, 753
# fri26.tsp best runs: 937

res_spm = bboptimize(make_perm_mapping_fitness(TSP, SPM);
    SearchRange = (0.0, 10.0), NumDimensions = size(TSP),
    PopulationSize = 5000,
    MaxTime = 30.0)
# dantzig42.tsp Best runs: 753

res_lrv = bboptimize(make_perm_mapping_fitness(TSP, LRV);
    SearchRange = (0.0, 10.0), NumDimensions = size(TSP),
    PopulationSize = 5000,
    MaxTime = 30.0)
# dantzig42.tsp Best runs: 881

res_spm = bboptimize(make_perm_mapping_fitness(TSP, SPM);
    SearchRange = (0.0, 10.0), NumDimensions = size(TSP),
    Method = :xnes,
    MaxTime = 30.0)
# dantzig42.tsp Best runs: 853

res_trc1000 = bboptimize(make_perm_mapping_fitness(TSP, TRC);
    SearchRange = searchrange(TRC), NumDimensions = size(TSP),
    PopulationSize = 1000,
    MaxTime = 30.0)
# dantzig42.tsp Best runs: 699, 722, 699, 699
# fri26.tsp best runs: 937

res_trc100 = bboptimize(make_perm_mapping_fitness(TSP, TRC);
    SearchRange = searchrange(TRC), NumDimensions = size(TSP),
    PopulationSize = 100,
    MaxTime = 30.0)
# dantzig42.tsp Best runs: 737, 778

# The TRC mapping is good on this problem since it "repairs" duplicates
# by turning them into the closest remaining city in the natural encoding.

# An interesting representation might be to encode the probability 
# to take the road from one city to another. This scales as O(n^2) but would 
# enable global learning since the same position in the genome means the same
# thing for all individuals. We only add the starting city as a rounded float
# and then use the highest probability of the already non-visited nodes
# from then on to create a unique TSP traversal.
struct EdgeTraversalProbabilitiesMapping <: PermutationMapping
    N::Int
    EdgeTraversalProbabilitiesMapping(N::Int) = new(N)
end
function searchspace(m::EdgeTraversalProbabilitiesMapping)
    ss = [(1.0, float(m.N)+0.99999999)]
    for i in 1:(m.N^2)
        push!(ss, (0.0, 1.0))
    end
    ss
end

apply(m::EdgeTraversalProbabilitiesMapping, v::Vector{Float64}) =
    apply!(m, v, Array{Int}(undef, m.N))

function apply!(m::EdgeTraversalProbabilitiesMapping, 
    v::Vector{Float64}, permutation::Vector{Int})

    i = 1
    permutation[i] = currentcity = floor(Int, v[1])
    visited = BitArray{1}(zeros(m.N))
    visited[currentcity] = 1
    while i < m.N
        i += 1
        startidx = 2 + (currentcity-1)*m.N
        endidx = startidx + m.N - 1
        perm = sortperm(view(v, startidx:endidx), rev=true)
        for candidatecity in perm
            if candidatecity != currentcity && !visited[candidatecity]
                permutation[i] = candidatecity
                visited[Int(candidatecity)] = 1
                currentcity = candidatecity
                break
            end
        end
    end
    permutation
end

const ETP = EdgeTraversalProbabilitiesMapping(size(TSP))

# Doesn't seem to work well, I guess the scaling makes it costly.
res_etp1000 = bboptimize(make_perm_mapping_fitness(TSP, ETP);
    SearchSpace = searchspace(ETP),
    PopulationSize = 1000,
    MaxTime = 60.0)

# We could try TRC with backups i.e. each also has a backup which is checked
# and only if that also is not available do we search for the closest remaining one.