"""
Base abstract class for black-box optimization algorithms.
"""
abstract type Optimizer end

"""
Optimizers derived from `SteppingOptimizer` implement classical iterative optimization scheme
`step!()` → `step!()` → …
"""
abstract type SteppingOptimizer <: Optimizer end
evaluator(so::SteppingOptimizer) = so.evaluator

"""
    step!(opt::SteppingOptimizer)

Do one iteration of the method.
"""
function step! end # FIXME avoid defining 0-arg function

"""
Base abstract class for optimizers that perform
`ask()` → ..eval fitness.. → `tell!()`
sequence at each step.
"""
abstract type AskTellOptimizer <: Optimizer end

"""
    ask(ato::AskTellOptimizer)

Ask for a new candidate solution to be generated, and a list of individuals
it should be ranked with.

The individuals are supplied as an array of tuples
with the individual and its index.

See also `tell!()`
"""
function ask end # FIXME avoid defining 0-arg function

# Our design is inspired by the object-oriented, ask-and-tell "optimizer API
# format" as proposed in:
#
#  Collette, Y., N. Hansen, G. Pujol, D. Salazar Aponte and
#  R. Le Riche (2010). On Object-Oriented Programming of Optimizers -
#  Examples in Scilab. In P. Breitkopf and R. F. Coelho, eds.:
#  Multidisciplinary Design Optimization in Computational Mechanics, Wiley,
#  pp. 527-565.
#  https://www.lri.fr/~hansen/collette2010Chap14.pdf
#
# but since Julia is not OO this is more reflected in certain patterns of how
# to specify and call optimizers. The basic ask-and-tell pattern is:
#
#   while !optimizer.stop
#     x = ask(optimizer)
#     y = f(x)
#     optimizer = tell(optimizer, x, y)
#   end
#
# after which the best solutions can be found by:
#
#   yopt, xopt = best(optimizer)
#
# We have extended this paradigm with the use of an archive that saves
# information on what we have learnt about the search space as well as the
# best solutions found. For most multi-objective optimization problems there
# is no single optimum. Instead there are many pareto optimal solutions.
# An archive collects information about the pareto optimal set or some
# approximation of it. Different archival strategies can be implemented.

"""
    tell!(ato::AskTellOptimizer, rankedCandidates)

Tell the optimizer about the ranking of candidates.
Returns the number of `rankedCandidates` that were inserted into the population,
because of the improved fitness.

See also `ask()`.
"""
function tell! end # FIXME avoid defining 0-arg function

"""
Base class for population-based optimization methods.
"""
abstract type PopulationOptimizer <: AskTellOptimizer end

population(popopt::PopulationOptimizer) = popopt.population
popsize(popopt::PopulationOptimizer) = popsize(population(popopt))

function setup!(o::SteppingOptimizer)
    # Do nothing, override if you need to setup prior to the optimization loop
end

function setup!(o::AskTellOptimizer, evaluator::Evaluator)
    # Do nothing, override if you need to setup prior to the optimization loop
end

function shutdown!(o::Optimizer)
    # Do nothing, override if you need to setup prior to the optimization loop
end

shutdown!(o::SteppingOptimizer) =
  shutdown!(evaluator(o)) # shutdown the evaluator

# The standard name function converts the type of the optimizer to a string
# and strips off trailing "Opt".
function name(o::Optimizer)
    s = string(typeof(o))
    if s[end-2:end] == "Opt"
        return s[1:end-3]
    else
        return s
    end
end

"""
    trace_state(io::IO, optimizer::Optimizer, mode::Symbol)

Trace the current optimization state to a given IO stream.
Called by `OptRunController` `trace_progress()`.

Override it for your optimizer to generate method-specific diagnostic traces.
"""
function trace_state(io::IO, optimizer::Optimizer, mode::Symbol) end
