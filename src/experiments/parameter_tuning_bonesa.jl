"""
A `BonesaTuner` has state for optimizing a certain algorithm on a set of
parameters described by a ParameterSet. The parameter vectors that are tuned
only take values in [0.0, 1.0] and it is up to the ParameterSet to translate
that to meaningful ranges.
"""
struct BonesaTuner
    numvectors::Int                     # Current num of param vectors in archive
    parameter_vectors::Matrix{Float64}  # Archive of param vectors, each column is one vector
    utilities::Matrix{Float64}          # Each column has Nu * Np utility values for one and the same parameter vector
    num_utility_values
end

# The predicted utility of a parameter vector x on problem is a weigthed, kernel
# smoothed average of the utilities of the vectors in the archive.
function predicted_utility(ba::BonesaArchive, x, problem)
end

function time_is_up(bt::BonesaTuner)
    time() - bt.start_time >= bt.max_time_seconds
end

# Start tuning parameters given a maximum time budget in hours.
function tune_parameters(bt::BonesaTuner, max_hours = 1.0)
    # Set up for measuring time left.
    bt.start_time = time()
    bt.max_time_seconds = max_hours * 60.0 * 60.0

    # Create a random set of vectors and then evaluate them on the problem set.
    create_random_parameter_vectors(bt, bt.M)
    bt.utilities = zeros(num_utilities(bt) * num_problems(bt.problem_set), bt.M)
    for i in 1:num_param_vectors(bt)
        if time_is_up(bt)
            # Cut the size of the archive since we did not have time to eval them all... :(
            bt.parameter_vectors = bt.parameter_vectors[:,1:(i-1)]
            bt.utilities = bt.utilities[:,1:(i-1)]
            break
        end
        bt.utilities[:, i] = evaluate_param_vector(bt, bt.parameter_vectors[:,i])
    end

    # Iteratively add to archive as we learn more about which param vectors are good.
    while( !time_is_up(bt) )
        # Update constants and temp calcs needed to calc pareto strenght below
    end

    # Select the terminal set (pareto front) of the param vectors.
end

# Calc predicted utilities, variance and support as well as pareto dominance
# and strengths for all param vectors in the archive.
function calc_temp_values_for_archive(bt::BonesaTuner)
    # Update magic constant c based on size of archive
    bt.c = calc_approximate_magic_constant_c(size_of_archive(bt), numdims(bt))

    # Calc the good vectors (saved as (boolean) indicator values)
    mean_actual_utilities = mean(bt.utilities, 2) # Mean per row => mean for each utility value
end

function support_of_vector(bt::BonesaTuner, pvector)
end

# For now we use random sampling but use lhs here instead, long-term.
function create_random_parameter_vectors(bt::BonesaTuner, num_vectors)
    bt.parameter_vectors = rand(num_params(bt), num_vectors)
end

function num_utility_values(bt::BonesaTuner)
    if bt.num_utility_values == false
        # Evaluate a param vector once to know how many utility values are returned.
        res = evaluate_param_vector_on_problem(bt, bt.problem_set[1], rand(num_params(bt)))
        bt.num_utility_values = length(res) + 1 # Add one for time
    else
        bt.num_utility_values
    end
end

function evaluate_param_vector(bt, pvector)
    parameters = map_to_parameter_values(bt.parameter_set, pvector)
    nu = num_utility_values(bt)
    np = num_problems(bt.problem_set)
    utilities = zeros(np*nu)

    for i in 1:np
        p = bt.problem_set[i]
        start = 1 + (i-1) * nu
        tic()
        utilities[ustart:(ustart+nu-2)] = evaluate_param_vector_on_problem(bt, p, parameters)
        utilities[ustart+nu-1] = toq()
    end

    return utilities
end

# Find an approximation to the BONESA magic constant c. It is the value
# that makes sum(w) == 50, where w is exp(c * normalized_distance(x, y)) for
# m uniformly, randomly sampled on the l-dimensional hypercube.
function calc_approximate_magic_constant_c(m, l, num_repeats = 5, tolerance = 1e-3)
    # Do a number of random approximations and then average
    mean([find_approximate_c_given_random_sample_of_points(m, l, tolerance) for i in 1:num_repeats])
end

function find_approximate_c_given_random_sample_of_points(m, l, tolerance = 1e-3)
    # Randomly sample m l-dimensional points
    points = rand(l, m)

    # To calc the distances we need the std dev per dimension
    std_devs = std(points, 2)

    # Now calc the distances between points. We only sample them if m is large.
    distances = zeros(Float64, m, m)
    for i in 1:m
        x = points[:,i]
        for j in i:m
            y = points[:,j]
            distances[i, j] = distances[i, j] = std_distance(x, y, std_devs)
        end
    end

    # Find a value of c by using bisection so that we are close to 50.0
    f = (c) -> average_sum_of_distances_for_c(c, distances) - 50.0
    bisection(f, -100.0, 100.0, tolerance)
end

# Find a value v that makes abs(f(v)) <= tolerance where a and b are guesses
# for v where b > a. This implementation assumes that we can find a negative
# (or positive) value of the function by stepping down (or up) => will not
# work for any function but throws an ArgumentError after many iterations of trying.
function bisection(f, a, b = 2*a, tolerance = 1e-5, max_iterations = 1e4)
    if f(a) > 0 && f(b) > 0
        a = step_value_until_changes_sign(f, a, -1.0, max_iterations, max_iterations)
    elseif f(a) < 0 && f(b) < 0
        b = step_value_until_changes_sign(f, b, 1.0, max_iterations, max_iterations)
    end

    p = (a + b)/2
    error = abs(f(p))
    iterations = 0
    while iterations < max_iterations && error > tolerance
        iterations += 1
        (f(a)*f(p) < 0) ? (b = p) : (a = p)
        p = (a + b)/2
        error = abs(f(p))
    end
    if iterations >= max_iterations
        throw(ArgumentError("Could not find a good approximation after $(iterations) iterations."))
    else
        return p
    end
end

function step_value_until_changes_sign(func, value::Float64, step_size::Float64 = 1.0, max_iterations = 1e3)
    wanted_sign = -1 * sign(func(value))
    iterations = 0
    while iterations < max_iterations && sign(func(value)) != wanted_sign
        iterations += 1
        value = value + step_size
        step_size *= 2
    end
    if iterations >= max_iterations
        throw(ArgumentError("Could not invert the sign after $(iterations) iterations."))
    else
        return value
    end
end

function average_sum_of_distances_for_c(c, distances)
    mean( exp(c * collect(triu(distances))) )
end

function std_distance(x, y, std_devs)
    sqrt(1/l * sum(((x .- y) ./ std_devs).^2))
end
