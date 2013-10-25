GlobalOptim.jl
==============

GlobalOptim is a (experimental, work-in-progress) global optimization framework for Julia (http://julialang.org/). It supports both multi- and single-objective optimization problems. The focus is on meta-heuristic algorithms (DE, PSO, CMA-ES etc) but currently only Differential Evolution (DE) is implemented. 

Eventually we hope to provide a JuMP interface to these global optimizers but since it is not clear if JuMP supports multiple objectives it's not clear if it makes sense.
