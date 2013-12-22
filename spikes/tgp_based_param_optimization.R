# We will use Treed Gaussian Processes to optimize the parameter settings
# of a Compass Search implemented with GSS.
library(tgp)

#####################################################################
# 1. Optimize single run for highest success ratio with fewest fevals
#####################################################################

# Set some constants.
n <- 5 # Start with a 5-D problem

# 1a. First define the parameters and their values
#  max_seconds = 0.1*n - 10*n
#  step_size   = 0.001 - deltabox
#  random_order = [0, 1] # bool, we round it to get truth value
#  freq_adapt_order = [0, 1] # bool, we round it to get truth value

# 1b. Latin Hypercube sampling to get initial points
library(lhs)
num_initial_points <- 20
sample_param_space <- function(n, k) {
  ranges <- matrix(c(
    10.0*n, 0.1*n,          # max_seconds
    (10.0 - -10.0), 0.001,  # step_size
    1, 0,                   # random_order
    1, 0                    # freq_adapt_order
    ), byrow=TRUE, nrow=4)
  delta <- ranges[,1] - ranges[,2]
  s <- improvedLHS(4, k, 10)
  x <- ranges[,2] + delta * s
  x[3,] <- round(x[3,])
  x[4,] <- round(x[4,])
  x
}

x = sample_param_space(n, num_initial_points)
