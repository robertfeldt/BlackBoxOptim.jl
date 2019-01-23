# Design check: Inserting initial candidate/solution as starting point for optimization
Date: 2018-09-08

## Possible designs
To set an initial solution we can either:

1. Pass it along as a parameter and all the init routines need to use the given one instead of sampling a (random) one.

2. Let them setup as now and then "insert" a candidate solution at a later point in time (but before starting the actual optimization).

3. Instead of as a parameter it could be given as a parameter to bbsetup/bboptimize. But most probably this should be handled in bbsetup/bboptimize and then added as a parameter.

3 can utilize either 1 or 2 as a way to enable this.

## Check of existing ways that optimizers are initialized

### GeneratingSetSearch, ResamplingMemeticSearcher, SimultaneousPerturbationSA2
 - Currently: Calls rand_individual(ss) in the constructor. 
 - For 1: Could take a starting point in instead. Maybe one has to give a starting point and the bbsetup method instead ensures there is one given or generates one by sampling in the search space.
 - For 2: Would be as easy as setting the x field later.

### RandomSearcher
 - Doesn't really make sense for this one. Probably do nothing. Throwing an exception or warning might be bad for compare_optimizers etc when one really wants to do the same thing for all optimizers.

### SeparableNESOpt, XNESOpt, DXNESOpt
 - Currently: Has an ini_x parameter to constructor. Starts from random point unless specified.

### DiffEvoOpt
 - Currently: Calls population(problem, opts) to create the population. 
 - So if opts includes the parameter for the InitialSolution or InitCandidate the population function needs to include it.

### BorgMOEA
 - Currently: Similar to DiffEvoOpt but calls population(problem, opts, nafitness(IndexedTupleFitness{N,F}), ntransient=1).
 - Solution should be same as for DiffEvoOpt but ensure all variants of population method are fixed.