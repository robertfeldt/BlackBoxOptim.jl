# Invoke tgp to do Treed Gaussian Process regression for prediction and design
# of parameter experiments on black box optimizers.
# Copyright (c) 2013-2014 Robert Feldt, robert.feldt@gmail.com
library(tgp)
library(rjson)
library(lhs)


#####################################################################
# 0. Define constants and functions
#####################################################################
num_tgp_runs <- 4                 # Number of runs when building tgp models
improvedlhs_num_repeats <- 100    # Number of points to consider for each added design point with improvedLHS
lhs_sample_size <- 500            # Number of candidate points sampled when selecting design points
num_selected_design_points <- 30  # Number of selected design points for prediction

sample_param_space <- function(ranges, num_samples) {
  delta <- ranges[,2] - ranges[,1]
  s <- improvedLHS(nrow(ranges), num_samples, improvedlhs_num_repeats)
  t(ranges[,1] + delta * s)
}

save_matrix_to_json_file <- function(matrix, filepath) {
  fileConn<-file(filepath)
  writeLines(toJSON(matrix), fileConn)
  close(fileConn)
}


#####################################################################
# 1. Parse parameters from command line
#####################################################################

args <- commandArgs(trailingOnly = TRUE)

# Example args: 6 "c(0.0,3.0,0.0,3.0,0.0,2.0,0.0,2.0,0.6,1.0,-2.0,0.0)" 7 13 in.csv out.json

num_params <- as.integer(args[1]);
ranges <- matrix(eval(parse(text = args[2])), byrow=TRUE, nrow=num_params);
num_new_runs <- as.integer(args[3]);
response_column <- as.integer(args[4]); # Response that should be minimized
csvfile <- args[5];
outfile <- args[6];


#####################################################################
# 2. Read input csv
#####################################################################
runs <- read.csv(csvfile, header = TRUE)


#####################################################################
# 3. Build tgp model if there are any existing runs, if not we
#    generate initial points via lhs. Save points to out file.
#####################################################################
if(nrow(runs) > 0) {

  num_runs <- nrow(runs)

  # Extract the X design matrix. Its columns must be the first num_params
  # columns of the runs data frame.
  X = as.matrix(runs[,1:num_params])

  # Extract the Z response vector.
  Z = as.matrix(runs[,response_column])

  # Build model for selecting new points
  model <- btgp(X=X, Z=Z, pred.n=FALSE, R=num_tgp_runs)

  # And select the best points from a LHS sample.
  Xcandidates = sample_param_space(ranges, lhs_sample_size);
  num_points <- max(num_selected_design_points, num_new_runs);
  XX <- tgp.design(num_points, Xcandidates, model);

  # Now predict in those points
  pmodel <- btgp(X=X, Z=Z, XX=XX, corr="exp", improv=TRUE,
    Ds2x=TRUE, R=5, verb=0, krige = FALSE)

  # Predicted values
  ZZ = m1p$ZZ.mean

  # Find the num_new_points points with minimum values and select those points.
  mins = sort(ZZ)
  indmins = which(ZZ %in% mins[1:num_new_runs])
  XXsel <- XX[indmins,]

} else {

  XXsel = sample_param_space(ranges, num_params+1);

}

save_matrix_to_json_file(XXsel, outfile);
