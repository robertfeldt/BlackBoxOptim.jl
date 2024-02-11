using SpecialFunctions
#Copyright (c) 2010, Xin-She Yang
#All rights reserved.
# Trasnlation to Julia language, Fredrik bagge Carlson
#
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are met:
#
#* Redistributions of source code must retain the above copyright notice, this
#  list of conditions and the following disclaimer.
#
#* Redistributions in binary form must reproduce the above copyright notice,
#  this list of conditions and the following disclaimer in the documentation
#  and/or other materials provided with the distribution
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
#FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
#CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
#OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# """
#     cuckoo_search(f,X0;
# `lb = -Inf`,
# `ub = Inf`,
# `n = Number of nests (or different solutions)
# `pa = 0.25` Discovery rate of alien eggs/solutions
# `tol = 1.0e-5`,
# `max_iter = 1e5`,
# `timeout = Inf`)
#
# Change this if you want to get better results
# Based on implementation by
# @inproceedings{yang2009cuckoo,
#   title={Cuckoo search via L{\'e}vy flights},
#   author={Yang, Xin-She and Deb, Suash},
#   booktitle={Nature \& Biologically Inspired Computing, 2009. NaBIC 2009. World Congress on},
#   pages={210--214},
#   year={2009},
#   organization={IEEE}
# }
#
# """
function cuckoo_search(f, X0;
    lb=-Inf,
    ub = Inf,
    n = 25,
    pa = 0.25,
    tol = 1.0e-5,
    max_iter = 1e5,
    timeout = Inf)

    nd = length(X0) # n dims

    if lb == -Inf
        lb = X0-0.99999*abs.(X0)
    end
    if ub == Inf
        ub = X0+0.99999*abs.(X0)
    end

    # Random initial solutions
    nest = [copy(X0) for _ = 1:n] # n = num_nests
    nest[1] = X0
    for i=2:n
        nest[i] = @. lb + (ub-lb)*rand()
    end
    # Get the current best
    fitness = fill(Inf,n)
    fmin, bestnest, nest, fitness = update_nests!(f,nest,nest,fitness)
    new_nest = deepcopy(nest)
    N_iter=0
    ## Starting iterations
    t0 = time()
    while fmin > tol && N_iter < max_iter
        # Generate new solutions (but keep the current best)
        move_cuckoos!(nest,new_nest,bestnest,lb,ub)
        fnew, best, nest, fitness = update_nests!(f,nest,new_nest,fitness)
        # Update the counter
        N_iter += n
        if fnew < fmin
            fmin = fnew
            bestnest = best
        end
        # Discovery and randomization
        new_nest = empty_nests(nest,lb,ub,pa)
        # Evaluate this set of solutions
        fnew, best, nest, fitness = update_nests!(f,nest,new_nest,fitness)
        # Update the counter again
        N_iter += n
        # Find the best objective so far
        if fnew < fmin
            fmin = fnew
            bestnest = best
        end
        if time()-t0 > timeout
            @info("Cuckoo search: timeout $(timeout)s reached ($(time()-t0)s)")
            break
        end
    end ## End of iterations
    ## Post-optimization processing
    ## Display all the nests
    bestnest, fmin

end
## --------------- All subfunctions are list below ------------------
## Get cuckoos by ramdom walk
function move_cuckoos!(nest,new_nest,best,lb,ub)
    # Levy flights
    n = length(nest)
    # Levy exponent and coefficient
    # For details, see equation (2.21), Page 16 (chapter 2) of the book
    # X. S. Yang, Nature-Inspired Metaheuristic Algorithms, 2nd Edition, Luniver Press, (2010).
    β = 3/2
    βi = 1/β
    σ = (gamma(1+β)*sin(pi*β/2)/(gamma((1+β)/2)*β*2^((β-1)/2)))^(1/β)
    u,v,step = similar(best),similar(best),similar(best)
    for j = 1:n
        s = nest[j]
        # This is a simple way of implementing Levy flights
        # For standard random walks, use step=1;
        ## Levy flights by Mantegna’s algorithm
        u .= σ.*randn.()
        v .= randn.()
        @. step = u/abs(v)^βi
        # In the next equation, the difference factor (s-best) means that
        # when the solution is the best solution, it remains unchanged.
        @. step = 0.01*step.*(s-best)
        # Here the factor 0.01 comes from the fact that L/100 should the typical
        # step size of walks/flights where L is the typical lenghtscale;
        # otherwise, Levy flights may become too aggresive/efficient,
        # which makes new solutions (even) jump out side of the design domain
        # (and thus wasting evaluations).
        # Now the actual random walks or flights
        s += step .* randn.()
        # Apply simple bounds/limits
        new_nest[j] = clamp(s,lb,ub)
    end
end
## Find the current best nest
function update_nests!(f,nest,new_nest,fitness)
    # Evaluating all new solutions
    for j = 1:length(nest)
        fnew = f(new_nest[j])
        if fnew <= fitness[j]
            fitness[j] = fnew
            nest[j] = new_nest[j]
        end
    end
    # Find the current best
    fmin,K = findmin(fitness)
    best = nest[K]
    fmin,best,nest,fitness
end

## Replace some nests by constructing new solutions/nests
function empty_nests(nest,lb,ub,pa)
    # A fraction of worse nests are discovered with a probability pa
    n = length(nest)
    # Discovered or not -- a status vector
    K = rand(size(nest)) .> pa
    # In the real world, if a cuckoo’s egg is very similar to a host’s eggs, then
    # this cuckoo’s egg is less likely to be discovered, thus the fitness should
    # be related to the difference in solutions. Therefore, it is a good idea
    # to do a random walk in a biased way with some random step sizes.
    ## New solution by biased/selective random walks
    stepsize = rand().*(nest[randperm(n)] .- nest[randperm(n)])
    new_nest = nest .+ stepsize .* K
    for j = 1:length(nest)
        new_nest[j] = clamp(new_nest[j],lb,ub)
    end
    return new_nest
end

function Base.clamp(x::AbstractVector,l::AbstractVector,u::AbstractVector)
    map(x,l,u) do x,l,u
        clamp(x,l,u)
    end
end

using BenchmarkTools, StaticArrays
# 78.619 ms (861606 allocations: 61.70 MiB)
# 6.321 ms (146181 allocations: 12.43 MiB)
rosenbrock(x) =  (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
cuckoo_search(rosenbrock, [10.,10], tol=1e-9)
@btime cuckoo_search(rosenbrock, [10.,10.], lb=[-20,-20], ub=[20,20], tol=1e-9)
@btime cuckoo_search(rosenbrock, @SVector([10.,10.]), lb=@SVector([-20,-20]), ub=@SVector([20,20]), tol=1e-9)


@profiler cuckoo_search(rosenbrock, @SVector([10.,10.]), lb=@SVector([-20,-20]), ub=@SVector([20,20]), tol=1e-9)
