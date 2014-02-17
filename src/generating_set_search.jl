# Generating Set Search as described in Kolda2003:
#  Kolda, Tamara G., Robert Michael Lewis, and Virginia Torczon. "Optimization 
#  by direct search: New perspectives on some classical and modern methods." 
#  SIAM review 45.3 (2003): 385-482.
#

# GSS is a type of DirectSearch
abstract DirectSearcher <: Optimizer

# A direction generator generates the search directions to use at each step of 
# a GSS search.
abstract DirectionGenerator

type ConstantDirectionGen <: DirectionGenerator
  directions::Array{Float64, 2}

  ConstantDirectionGen(directions) = begin
    new(directions)
  end
end

function directions_for_k(cg::ConstantDirectionGen, k)
  cg.directions # For a ConstantDirectionGen it is always the same regardless of k...
end

# We can easily do a compass search with GSS by generating directions 
# individually (+ and -) for each coordinate.
compass_search_directions(n) = ConstantDirectionGen([eye(n,n) -eye(n, n)])

GSSDefaultParameters = {
  :DeltaTolerance => 1e-10,       # GSS has converged if the StepSize drops below this tolerance level
  :StepSizeFactor => 0.50,        # Factor times the minimum search space diameter to give the initial StepSize
  :RandomDirectionOrder => true,  # Randomly shuffle the order in which the directions are used for each step
}

type GeneratingSetSearcher <: DirectSearcher
  parameters::Parameters
  direction_gen::DirectionGenerator
  search_space::SearchSpace
  n::Int64
  k::Int64
  step_size::Float64
  x::Array{Float64, 2}
  xfitness::Float64

  GeneratingSetSearcher(parameters) = begin
    params = Parameters(parameters, GSSDefaultParameters)
    n = numdims(params[:Evaluator])
    ss = search_space(params[:Evaluator])
    dgen = get(params, :DirectionGenerator, compass_search_directions(n))
    step_size = params[:StepSizeFactor] * minimum(diameters(ss))
    x = rand_individual(ss)
    new(params, dgen, ss, n, 0, step_size, 
      x, evaluate(params[:Evaluator], x))
  end
end

has_ask_tell_interface(gss::GeneratingSetSearcher) = false

function has_converged(gss::GeneratingSetSearcher)
  gss.step_size < gss.parameters[:DeltaTolerance]
end

function step(gss::GeneratingSetSearcher)
  if has_converged(gss)
    # Restart from a random point
    gss.x = rand_individual(gss.search_space)
    gss.xfitness = evaluate(gss.parameters[:Evaluator], gss.x)
    gss.step_size = gss.parameters[:StepSizeFactor] * minimum(diameters(gss.search_space))
  end

  # Get the directions for this iteration
  gss.k += 1
  directions = directions_for_k(gss.direction_gen, gss.k)

  # Set up order vector from which we will take the directions after possibly shuffling it
  order = collect(1:size(directions, 2))
  if gss.parameters[:RandomDirectionOrder] == true
    shuffle!(order)
  end

  # Check all directions to find a better point; default is that no one is found.
  found_better = false
  candidate = zeros(gss.n, 1)

  # Loop over directions until we find an improvement (or there are no more directions to check).
  for(direction in order)

    candidate = gss.x + gss.step_size .* directions[:, direction]
    rand_bound_from_target!(candidate, gss.x, gss.search_space)

    if is_better(gss.parameters[:Evaluator], candidate, gss.xfitness)
      found_better = true
      break
    end

  end

  if found_better
    gss.x = candidate
    gss.xfitness = last_fitness(gss.parameters[:Evaluator])
  else
    gss.step_size = gss.step_size / 2
  end
end