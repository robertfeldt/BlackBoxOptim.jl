BlackBoxOptim.jl
==============

[![Build Status](https://travis-ci.com/robertfeldt/BlackBoxOptim.jl.svg?branch=master)](https://travis-ci.com/robertfeldt/BlackBoxOptim.jl)

`BlackBoxOptim` is a global optimization package for Julia (http://julialang.org/). It supports both multi- and single-objective optimization problems and is focused on (meta-)heuristic/stochastic algorithms (DE, NES etc) that do NOT require the function being optimized to be differentiable. This is in contrast to more traditional, deterministic algorithms that are often based on gradients/differentiability. It also supports parallel evaluation to speed up optimization for functions that are slow to evaluate.

# Installation
```julia
using Pkg; Pkg.add("BlackBoxOptim")
```
or latest master directly from github:
```julia
using Pkg; Pkg.clone("https://github.com/robertfeldt/BlackBoxOptim.jl")
```
from a Julia repl.

# Usage

To show how the `BlackBoxOptim` package can be used, let's implement the Rosenbrock function, a classic problem in numerical optimization. We'll assume that you have already installed `BlackBoxOptim` as described above.

First, we'll load `BlackBoxOptim` and define the Rosenbrock function (in 2 dimensions):
```julia
using BlackBoxOptim

function rosenbrock2d(x)
  return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
end
```
We can now call the `bboptimize()` function, specifying the function to be optimized (here: `rosenbrock2d()`) and the range of values allowed for each of the dimensions of the input:
```julia
res = bboptimize(rosenbrock2d; SearchRange = (-5.0, 5.0), NumDimensions = 2)
```
We get back an optimization result object that we can query to, for example, get the best fitness and solution candidate:
```julia
best_fitness(res) < 0.001        # Fitness is typically very close to zero here...
length(best_candidate(res)) == 2 # We get back a Float64 vector of dimension 2
```
`BlackBoxOptim` will default to using an adaptive differential evolution optimizer in this case and use it to try to locate a solution where both elements can be Floats in the range -5.0:5.0. If you wanted a different range of allowed values for the second dimension of the solution you can specify that with a range of allowed values. In this case you do not need to specify the number of dimensions since that is implicit from the number of ranges supplied:
```julia
bboptimize(rosenbrock2d; SearchRange = [(-5.0, 5.0), (-2.0, 2.0)])
```
If you want to use a different optimizer that can be specified with the `Method` keyword. For example, using the standard differential evolution optimizer `DE/rand/1/bin`:
```julia
bboptimize(rosenbrock2d; SearchRange = (-5.0, 5.0), NumDimensions = 2, Method = :de_rand_1_bin)
```
You can (this currently requires the master branch so `] add BlackBoxOptim#master`) give a starting (initial candidate) point for the search when calling `bboptimize` but beware
that very little checking is done on it so be sure to provide a candidate of the right length and 
inside the search space:
```julia
x0 = [1.0, 1.0] # starting point (aka initial candidate)
res = bboptimize(rosenbrock2d, x0; SearchRange = (-5.0, 5.0), NumDimensions = 2, MaxTime = 0.1)
isapprox(best_fitness(res), 0.0)
```

Note that the `rosenbrock2d()` function is quite easy to optimize. Even a random search will come close if we give it more time:
```julia
bboptimize(rosenbrock2d; SearchRange = (-5.0, 5.0), NumDimensions = 2, Method = :random_search, MaxTime = 10.0)
```
But if we optimize the same rosenbrock function in, say, 30 dimensions that will be very hard for a random searcher while sNES or DE can find good solutions if we give them some time. We can compare optimizers using the `compare_optimizers()` function:
```julia
function rosenbrock(x)
  return( sum( 100*( x[2:end] .- x[1:end-1].^2 ).^2 .+ ( x[1:end-1] .- 1 ).^2 ) )
end

res = compare_optimizers(rosenbrock; SearchRange = (-5.0, 5.0), NumDimensions = 30, MaxTime = 3.0);
```
You can find more examples of using `BlackBoxOptim` in [the examples directory](examples).

# Multi-objective optimization

Multi-objective evaluation is supported by the BorgMOEA algorithm. Your fitness function should return a tuple of the objective values and you should indicate the fitness scheme to be (typically) Pareto fitness and specify the number of objectives. Otherwise the use is similar, here we optimize the Schaffer1 function:
```julia
fitness_schaffer1(x) = (sum(abs2, x), sum(abs2, x .- 2.0))
res = bboptimize(fitness_schaffer1; Method=:borg_moea,
            FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true),
            SearchRange=(-10.0, 10.0), NumDimensions=3, ϵ=0.05,
            MaxSteps=50000, TraceInterval=1.0, TraceMode=:verbose);
```
`pareto_frontier(res)` would give a vector of all Pareto-optimal solutions and corresponding fitness values.
If we simply want to get one individual with the best aggregated fitness:
```julia
bs = best_candidate(res)
bf = best_fitness(res)
```
By default, the aggregated fitness is the sum of the individual objective values, but this could be changed when declaring the fitness scheme, e.g.
the weighted sum with weights `(0.3, 0.7)`:
```julia
weightedfitness(f) = f[1]*0.3 + f[2]*0.7

...
    FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true, aggregator=weightedfitness)
...
```
Of course, once the Pareto set (`pareto_frontier(res)`) is found, one
can apply different criteria to filter the solutions.
For example, to find the solution with the minimal first objective:
```julia
pf = pareto_frontier(res)
best_obj1, idx_obj1 = findmin(map(elm -> fitness(elm)[1], pf))
bo1_solution = params(pf[idx_obj1]) # get the solution candidate itself...
```
or to use different weighted sums:
```julia
weighedfitness(f, w) = f[1]*w + f[2]*(1.0-w)
weight = 0.4 # Weight on first objective, so second objective will have weight 1-0.4=0.6
best_wfitness, idx = findmin(map(elm -> weighedfitness(fitness(elm), weight), pf))
bsw = params(pf[idx])
```

# Configurable Options

The section above described the basic API for the `BlackBoxOptim` package. There is a large number of different optimization algorithms that you can select with the `Method` keyword (`adaptive_de_rand_1_bin`, `adaptive_de_rand_1_bin_radiuslimited`, `separable_nes`, `xnes`, `de_rand_1_bin`, `de_rand_2_bin`, `de_rand_1_bin_radiuslimited`, `de_rand_2_bin_radiuslimited`, `random_search`, `generating_set_search`, `probabilistic_descent`, `borg_moea`).

In addition to the `Method` parameter, there are many other parameters you can change. Some key ones are:

* `MaxTime`: For how long can the optimization run? Defaults to `false` which means that number of iterations is the given budget, rather than time.
* `MaxFuncEvals`: How many evaluations that are allowed of the function being optimized.
* `TraceMode`: How optimization progress should be displayed (`:silent`, `:compact`, `:verbose`). Defaults to `:compact` that outputs current number of fitness evaluations and best value each `TraceInterval` seconds.
* `PopulationSize`: How large is the initial population for population-based optimizers? Defaults to `50`.
* `TargetFitness`. Allows to specify the value of the best fitness for a given problem. The algorithm stops as soon as the distance between the current `best_fitness()` and `TargetFitness` is less than `FitnessTolerance`.
This list is not complete though, please refer to the `examples` and `tests` directories for additional examples.

# State of the Library

## Existing Optimizers

* Natural Evolution Strategies:
  - Separable NES: `separable_nes`
  - Exponential NES: `xnes`
  - Distance-weighted Exponential NES: `dxnes`
* Differential Evolution optimizers, 5 different:
  - Adaptive DE/rand/1/bin: `adaptive_de_rand_1_bin`
  - Adaptive DE/rand/1/bin with radius limited sampling: `adaptive_de_rand_1_bin_radiuslimited`
  - DE/rand/1/bin: `de_rand_1_bin`
  - DE/rand/1/bin with radius limited sampling (a type of trivial geography): `de_rand_1_bin_radiuslimited`
  - DE/rand/2/bin: `de_rand_2_bin`
  - DE/rand/2/bin with radius limited sampling (a type of trivial geography): `de_rand_2_bin_radiuslimited`
* Direct search:
  - Generating set search:
    - Compass/coordinate search: `generating_set_search`
    - Direct search through probabilistic descent: `probabilistic_descent`
* Resampling Memetic Searchers:
  - Resampling Memetic Search (RS): `resampling_memetic_search`
  - Resampling Inheritance Memetic Search (RIS): `resampling_inheritance_memetic_search`
* Stochastic Approximation:
  - Simultaneous Perturbation Stochastic Approximation (SPSA): `simultaneous_perturbation_stochastic_approximation`
* RandomSearch (to compare to): `random_search`

For multi-objective optimization only the [BorgMOEA](http://borgmoea.org/) (`borg_moea`) is supported but it is a good one. :)

# Multithreaded and Parallel Function Evaluation

**NB!** There are problems with the multi-threaded evaluation on Julia 1.6 and later. We will be investigating this and hope to fix in a future release. For now the related tests have been deactivated. Sorry for the inconvenience.

For some (slow) functions being optimized and if you have a multi-core CPU you can gain performance by using multithreaded or parallel evaluation. This typically requires an optimization algorithm that evaluates many candidate points in one batch. The NES family (`xnes`, `dxnes` etc), for example. See the file

[examples/rosenbrock_parallel.jl](examples/rosenbrock_parallel.jl)

for one example of this. On Julia 1.3 and later it is typically better to use multithreading see the file

[examples/multithreaded_optimization.jl](examples/multithreaded_optimization.jl)

for some examples.

# Guide to selecting an optimizer

In our experiments the radius limited DE's perform better than the classic `de_rand_1_bin DE` in almost all cases. And combining it with adaptive setting of the weights makes it even better. So for now `adaptive_de_rand_1_bin_radiuslimited()` is our recommended "goto" optimizer. However, the difference between the top performing DE's is slight.

Some NES variants (separable or dx) can sometimes beat the DE optimizers in the tests we have done. But xnes and dxnes are often very slow and while the separable NES isn't it is less robust. So we don't recommend it as a default, robust choice.

We maintain a [list of optimizers ranked by performance](examples/benchmarking/latest_toplist.csv) when tested on a large set of problems. From the list we can see that the adaptive differential evolution optimizers (`adaptive_de_rand_1_bin` and/or `adaptive_de_rand_1_bin_radiuslimited`) are typically on top when it comes to mean rank. The `generating_set_search` often gives best results (its `MedianLogTimesWorseFitness` is typically in the range `0.3`-`0.6`, which means its median fitness value is `10^0.3=2.0` to `10^0.6=4.0` times worse than the best fitness found on a problem) and is faster (ranked first on run time, typically) but it is not as robust as the DE optimizers and thus is ranked lower on mean rank (per problem).

Overall we recommend one of the DE optimizers as a generally good choice since their runtime is often good and they are robust and works well both for simpler, separable problems as well as for more complex ones. They also tend to scale better to high-dimensional settings. Of course, optimizer performance varies depending on the problem and the dimensionality so YMMV. 
