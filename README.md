GlobalOptim.jl
==============

GlobalOptim is a (experimental, work-in-progress) global optimization framework for Julia (http://julialang.org/). It supports both multi- and single-objective optimization problems and is focused on (meta-)heuristic/stochastic algorithms (DE, PSO, CMA-ES etc) rather than more traditional, deterministic algorithms (as available in the Optim.jl library).

Eventually we hope to provide a JuMP interface but since it is not clear if JuMP supports multiple objectives this is to be decided.

Usage
=====

TBD

State of the Library
====================

Existing Optimizers
-------------------

* Differential Evolution optimizers:
  - DE/rand/1/bin: de_rand_1_bin()
  - DE/rand/1/bin with radius limited sampling (a type of trivial geography): de_rand_1_bin_radiuslimited()
  - AdaptiveDE/rand/1/bin: de_rand_1_bin()
  - AdaptiveDE/rand/1/bin with radius limited sampling: adaptive_de_rand_1_bin_radiuslimited()

Planned Optimizers
------------------

* Particle Swarm, PSO
* CMA-ES
* Amalgam meta-optimizer (by Vrugt), which takes a set of (at least 2) other optimizers and switches between them dynamically during the search.

Utilities
---------
* Latin hypercube sampling for creating initial populations

Guide to selecting an optimizer
===============================

In our experiments the radius limited DE's perform better than the classic de_rand_1_bin DE in almost all cases. And combining it with adaptive setting of the weights makes it even better. So for now adaptive_de_rand_1_bin_radiuslimited() is our recommended "goto" optimizer. However, the difference between the top performing DE's is slight.

Once we have Amalgam implemented we believe that Amalgam(PSO, CMA-ES, DE/rand/1/bin/radiuslimited) will be a very powerful default choice. This remains to be evaluated though.