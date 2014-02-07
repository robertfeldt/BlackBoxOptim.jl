BlackBoxOptim.jl
==============

BlackBoxOptim is a (work-in-progress) global optimization framework for Julia (http://julialang.org/). It supports both multi- and single-objective optimization problems and is focused on (meta-)heuristic/stochastic algorithms (DE, NES, CMA-ES etc) that do NOT require the function being optimized to be differentiable. This is in contrast to more traditional, deterministic algorithms that are often based on gradients/differentiability.

Eventually we hope to provide a JuMP interface but since it is not clear if JuMP supports multiple objectives this is to be decided.

# Installation

Just install from the github by calling:

    Pkg.clone("https://github.com/robertfeldt/BlackBoxOptim.jl")

from a Julia repl.

# Usage

To show how the BlackBoxOptim package can be used, let's implement the Rosenbrock function, a classic problem in numerical optimization. We'll assume that you have already installed BlackBoxOptim as described above.

First, we'll load BlackBoxOptim and define the Rosenbrock function (in 2 dimensions):

    using BlackBoxOptim

    function rosenbrock2d(x)
      return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
    end

We can now call the `bboptimize` function, specifying the function to be optimized (here: rosenbrock2d) and the range of values allowed for each of the dimensions of the input:

    bboptimize(rosenbrock2d, (-5.0, 5.0); dimensions = 2)

BlackBoxOptim will default to using an adaptive differential evolution optimizer in this case and use it to try to locate a solution where both elements can be Floats in the range -5.0:5.0. If you wanted a different range of allowed values for the second dimension of the solution you can specify that with a range of allowed values. In this case you do not need to specify the number of dimensions since that is implicit from the number of ranges supplied:

    bboptimize(rosenbrock2d, [(-5.0, 5.0), (-2.0, 2.0)])

If you want to use a different optimizer that can be specified with the `method` keyword. For example, using the standard differential evolution optimizer DE/rand/1/bin:

    bboptimize(rosenbrock2d, (-5.0, 5.0); dimensions = 2, method = :de_rand_1_bin)

Note that the rosenbrock2d function is quite easy to optimize. Even a random search will come close if we give it more time:

    bboptimize(rosenbrock2d, (-5.0, 5.0); dimensions = 2, method = :random_search, max_time = 10.0)

But if we optimize the same rosenbrock function in, say, 30 dimensions that will be very hard for a random searcher while sNES or DE can find good solutions if we give them some time. We can compare optimizers using the `compare_optimizers` function:

    function rosenbrock(x)
      return( sum( 100*( x[2:end] - x[1:end-1].^2 ).^2 + ( x[1:end-1] - 1 ).^2 ) )
    end

    res = compare_optimizers(rosenbrock, (-5.0, 5.0); dimensions = 30, max_time = 5.0);

# Configurable Options

The section above described the basic API for the BlackBoxOptim package. We employed several different optimization algorithms using the `method` keyword, which can take on any of the following values:

* `separable_nes`
* `xnes`
* `de_rand_1_bin`
* `de_rand_2_bin`
* `de_rand_1_bin_radiuslimited`
* `adaptive_de_rand_1_bin`
* `adaptive_de_rand_1_bin_radiuslimited`
* `random_search`

In addition to the `method` keyword, you can alter the behavior of the Optim package by using other keywords:

* `iterations`: How many iterations (function evaluations) will be run before the algorithm gives up? Defaults to `10_000`.
* `max_time`: For how long can the optimization run? Defaults to false which means that number of iterations is the given budget, rather than time.
* `store_trace`: Should a trace of the optimization be stored? Defaults to `false`.
* `show_trace`: Should a trace of the optimization be shown on `STDOUT`? Defaults to `false`.
* `population_size`: How large is the initial population for population-based optimizers? Defaults to `50`.
* `method_options`: Detailed control over the optimization method, through a Dict mapping named options to their values. Defaults to {} which means that the default options are used for each method.

Most optimizers have specific options that can be specified in the `method_options` dict. Further details TBD.

# State of the Library

## Existing Optimizers

* Natural Evolution Strategies:
  - Separable NES: `separable_nes()`
  - Exponential NES: `xnes()`
* Differential Evolution optimizers, 5 different:
  - DE/rand/1/bin: `de_rand_1_bin()`
  - DE/rand/1/bin with radius limited sampling (a type of trivial geography): `de_rand_1_bin_radiuslimited()`
  - DE/rand/2/bin: `de_rand_2_bin()`
  - Adaptive DE/rand/1/bin: `de_rand_1_bin()`
  - Adaptive DE/rand/1/bin with radius limited sampling: `adaptive_de_rand_1_bin_radiuslimited()`
* RandomSearch (to compare to): `random_search()`

## Planned Optimizers

* HillClimbing
* CMA-ES
* Accelerated Coordinate Descent
* Amalgam meta-optimizer (by Vrugt), which takes a set of (at least 2) other optimizers and switches between them dynamically during the search.

## Utilities
* Latin hypercube sampling for creating initial populations
* Frequency Adaptation (of methods/strategies/coordinates)

## Planned Utilities
* Comparison of optimizers based on (Moré's') "data/performance profiles".
* Running BBOB/COCO comparisons of optimizers

## Problems

* Rosenbrock
* Sphere
* Schwefel2.21
* Schwefel2.22
* Schwefel1.2
* Step
* Rastrigin
* Ackley
* Griewank
* Ellipsoid

...and more, see https://github.com/robertfeldt/BlackBoxOptim.jl/blob/master/src/problems/single_objective.jl
for details.

## Planned Problems
* BBOB/COCO (Black-Box Optimization Benchmark / COmparing Continuous Optimizers) problems:
  - Separable (Unimodal):
    1. Sphere
  - Low or moderate condition (Unimodal):
    8. Rosenbrock, original
  - High condition (Unimodal):
    10. Ellipsoid with monotone x-transformation, condition 1e6
  - Multi-modal:
    15. Rastrigin with both x-transformations, condition 10

# Guide to selecting an optimizer

In our experiments the radius limited DE's perform better than the classic de_rand_1_bin DE in almost all cases. And combining it with adaptive setting of the weights makes it even better. So for now adaptive_de_rand_1_bin_radiuslimited() is our recommended "goto" optimizer. However, the difference between the top performing DE's is slight.

The separable NES often beats all of the DE optimizers in the tests we have done. But it is about 2-3 times slower per iteration so not really a fair comparison. It seems it can still hold up even if we normalize for time rather than number of executions but since it is not as good for non-separable problems it is not our default for now. XNES can sometimes beat sNES but scales very badly so is not a good default choice.

Once we have Amalgam implemented we believe that Amalgam(PSO, CMA-ES, DE/rand/1/bin/radiuslimited) will be a very powerful default choice. This remains to be evaluated though.