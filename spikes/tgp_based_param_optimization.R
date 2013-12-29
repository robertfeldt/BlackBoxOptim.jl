# We will use Treed Gaussian Processes to optimize the parameter settings
# of a Compass Search implemented with GSS.
library(tgp)

#####################################################################
# 1. Optimize single run for highest success ratio with fewest fevals
#####################################################################

# Set some constants.
n <- 2 # Start with a simple 2-D cigar problem on [-10, 10]^n => diameter = 20
diameter <- 10 - (-10)
problem <- "cigar"

# 1a. First define the parameters and their values
#  max_seconds      = 10^[-2, 1]*n
#  step_size        = 10^[-2, 0]*sample_space_diameter
#  random_order     = [0, 1]        # bool, we round it later
#  freq_adapt_order = [0, 1]        # bool, we round it later
#  delta_tol        = 10^[-13, -4]

# 1b. Latin Hypercube sampling to get initial points
library(lhs)
num_params = 5

sample_param_space <- function(ranges, num_samples) {
  delta <- ranges[,2] - ranges[,1]
  s <- improvedLHS(num_samples, nrow(ranges), 100)
  ranges[,1] + delta * s
}

experiment_ranges <- matrix(c(
    -1, 1,          # log(max_seconds)
    -2, 1,          # log(step_size)
     0,  1,         # random_order
     0,  1,         # freq_adapt_order
    -13, -2         # log(delta_tol)
), byrow=TRUE, nrow=num_params);

x = sample_param_space(experiment_ranges, num_params+1)

#           [,1]       [,2]         [,3]
#[1,]  0.7146809 -0.4500064   0.03894958
#[2,]  0.2793796 -0.1929759   0.80923724
#[3,]  0.3368907  0.1240581   0.68435067
#[4,]  0.1799453  0.8234069   0.13939124
#[5,] -8.3559134 -6.8995563 -10.68603265

gen_julia_call <- function(params, problem) {
  cat(paste("xb, fb, fevals, reason = generating_set_search(", problem, ", ", n, "; ",
    "max_seconds = ", n*10^params[1], 
    ", delta_tol = ", 10^params[5], 
    ", step_size = ", diameter*10^params[2], 
    ", x = false, random_order = ", ifelse(round(params[3])==1.0, "true", "false"),
    ", freq_adapt_order = ", ifelse(round(params[4])==1.0, "true", "false"),
    ")", sep=""))
}

xbs = matrix(rep(0.0, n*3), nrow=3)
y = rep(0.0, 3)
fevals = rep(0.0, 3)

# Run 1:
gen_julia_call(x[,1], problem)
# generating_set_search(cigar, 2; max_seconds = 0.110843199260581, delta_tol = 1.24047622585845e-12, step_size = 3.84433676126225, x = false, random_order = true, freq_adapt_order = false)
xbs[1,] = c(8.21044e-14, 3.99243e-13)
y[1] = 1.593949549092762e-19
fevals[1] = 251
reason1 = "Step size smaller than delta_tol"

# Run 2:
gen_julia_call(x[,2], problem)
# generating_set_search(cigar, 2; max_seconds = 0.0379017703267553, delta_tol = 5.53792388484792e-05, step_size = 2.37168525931216, x = false, random_order = false, freq_adapt_order = true)
xbs[2,] = c(-2.62615e-5, -2.92947e-5)
y[2] = 0.0008581784581966408
fevals[2] = 106
reason2 = "Step size smaller than delta_tol"

# Run 3:
gen_julia_call(x[,3], problem)
# generating_set_search(cigar, 2; max_seconds = 0.586566493452953, delta_tol = 1.2024516664737e-06, step_size = 15.8482815318414, x = false, random_order = false, freq_adapt_order = false)
xbs[3,] = c(6.82709e-7, -3.65295e-7)
y[3] = 1.334409736334325e-7
fevals[3] = 172
reason3 = "Step size smaller than delta_tol"

# Now build the first tgp model!
XX = t(improvedLHS(num_params, 100, 100))
m1_btgp <- btgp(X=t(x), XX=XX, Z = y)

# GP paket att testa: tgp, dynatree, GPfit, mlegp