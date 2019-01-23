# This is based on the COCO (COmparing Continuous Optimizers) source code
# and tries to stay close to its methods/functions in order to allow the same
# type of comparisons in Julia.
module COCO

# BBOB test functions are optimization problems
abstract type BBOBFunction <: OptimizationProblem end

# Number of dimensions.
ndims(f::BBOBFunction) = f.ndims

instantiated(f::BBOBFunction) = ndims(f) != Nothing

# Returns the best function value of this instance of the function.
getfopt(f::BBOBFunction) = instantiated(f) ? f.fopt : Exception("This function class has not been instantiated yet.")

# Set/get the argument of the optimum of the function.
setxopt(f::BBOBFunction, xopt) = f.xopt = xopt
getxopt(f::BBOBFunction) = f.xopt

# Returns the objective function value.
evaluate(f::BBOBFunction, x) = evalfull(f, x)[0]

# Returns a float penalty for being outside of boundaries [-5, 5]
function defaultboundaryhandling(x, fac = 1.0)
    xoutside = maximum(vcat(zeros(size(x)), abs(x) - 5), 1) .* sign(x)
    return fac * sum(xoutside.^2)
end

# Most functions do not have a penalty outside the boundaries though so return
# a zero penalty by default.
boundaryhandling(f::BBOBFunction, x) = 0

function evalfull(f::BBOBFunction, x)
    fadd = getfopt(f)
    curshape, dim = size(x)
    # it is assumed x are row vectors

    if f.lastshape != curshape
      initwithsize(curshape, dim)
    end

    # BOUNDARY HANDLING
    fadd += boundaryhandling(f, x)

    # TRANSFORMATION IN SEARCH SPACE
    x = x - arrxopt(f) # cannot be replaced with x -= arrxopt! RF: Why???

    # COMPUTATION core
    ftrue = compute_core(f, x)
    fval = noise(f.noisefunc, ftrue)

    # FINALIZE
    ftrue += fadd
    fval += fadd
    return (fval, ftrue)
end

# All BBOBFunction functions must implement evalfull, which return noisy and
# noise-free value, the latter for recording purpose.
# evalfull(f::BBOBFunction, x) = ...

# Boundary handling
# boundaryhandling(f::BBOBFunction, x) = ...

# We depart from the BBOB type/class hierarchy they use in Python since
# we can use Julia's type system to model the noise functions in a better
# way than how it is done in Python. Thus a noise function need not be a
# BBOBFunction.
abstract type BBOBNoiseFunction end

# Different types depending on the noise function used.
struct BBOBNfreeFunction <: BBOBNoiseFunction
end
noise(f::BBOBNfreeFunction, ftrue) = copy(ftrue) # no noise added

struct BBOBGaussFunction <: BBOBNoiseFunction
    gaussbeta::Float64
end
noise(f::BBOBGaussFunction, ftrue) = fGauss(ftrue, f.gaussbeta)

struct BBOBUniformFunction <: BBOBNoiseFunction
    unifalphafac::Float64
    unifbeta::Float64
end
noise(f::BBOBUniformFunction, ftrue) =
    fUniform(ftrue, f.unifalphafac * (0.49 + 1. / f.dim), f.unifbeta)

struct BBOBCauchyFunction <: BBOBNoiseFunction
    cauchyalpha::Float64
    cauchyp::Float64
end
noise(f::BBOBCauchyFunction, ftrue) = fCauchy(ftrue, f.cauchyalpha, f.cauchyp)

abstract type FSphere{NoiseFunc} <: BBOBFunction end
compute_core(f::FSphere, x) = sum(x.^2)

# Sphere without noise
struct F1 <: FSphere{BBOBNfreeFunction}
    funId::Int
    noisefunc::BBOBNoiseFunction
    F1() = new(1, BBOBNfreeFunction())
end
boundaryhandling(f::F1, x) = 0

# Sphere with Gaussian noise
struct F101 <: FSphere{BBOBGaussFunction}
    funId::Int
    noisefunc::BBOBNoiseFunction
    F101 = new(101, BBOBGaussFunction(0.01))
end

# Run an optimizer like in the COCO (COmparing Continuous Optimizers) sw,
# i.e. with a given number of dimensions, specific target value and max number
# of function evaluations.
# COCO method: MY_OPTIMIZER in MY_OPTIMIZER.m
# COCO differences:
#  - We send in an optimizer instead of changing the code each time we want to evaluate an optimizer.
#function coco_run_optimizer(optimizer, problem, numDimensions,
#  funcTarget, maxFunevals; populationSize = 200)
#  maxfunevals = min(1e8 * numDimensions, maxFunevals)
#  popsize = min(maxFunevals, populationSize)
#  fbest = inf
#  for iter in 1:ceil(maxfunevals/popsize)
#
#    xpop = 10 * rand(numDimensions, popsize) - 5;      % new solutions
#    [fvalues, idx] = sort(feval(FUN, xpop)); % evaluate
#    if fbest > fvalues(1)                    % keep best
#      fbest = fvalues(1);
#      xbest = xpop(:,1);
#    end
#    if feval(FUN, 'fbest') < ftarget         % COCO-task achieved
#      break;                                 % (works also for noisy functions)
#    end
#  end
#end
#
#function coco
