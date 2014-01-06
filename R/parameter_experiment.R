# Invoke tgp to do Treed Gaussian Process regression for prediction and design
# of parameter experiments on black box optimizers.
# Copyright (c) 2013-2014 Robert Feldt, robert.feldt@gmail.com
library(tgp)
library(rjson)
library(lhs)
library(optparse, quietly = TRUE, warn.conflicts = FALSE, verbose = FALSE)


#####################################################################
# 0. Define constants and functions
#####################################################################
num_tgp_runs <- 4                 # Number of runs when building tgp models
improvedlhs_num_repeats <- 10     # Number of points to consider for each added design point with improvedLHS
lhs_sample_size <- 500            # Number of candidate points sampled when selecting design points
num_selected_design_points <- 50  # Number of selected design points for prediction

sample_param_space <- function(num_params, num_samples, repeats = improvedlhs_num_repeats) {
  #delta <- ranges[,2] - ranges[,1]
  s <- improvedLHS(num_samples, num_params, repeats)
  #t(ranges[,1] + delta * s)
  s
}

save_matrix_to_json_file <- function(matrix, filepath) {
  fileConn<-file(filepath)
  writeLines(toJSON(matrix), fileConn)
  close(fileConn)
}

save_list_to_json_file <- function(list, filepath) {
  fileConn<-file(filepath)
  writeLines(toJSON(list), fileConn)
  close(fileConn)
}


#####################################################################
# 1. Parse parameters from command line
#####################################################################

args <- commandArgs(trailingOnly = TRUE)

# Example 
# args <- c("6", "1", "15", "4", "cmsa_es_exp2.csv", "new_runs.json", "sa")
# args <- c("4", "1", "-11", "4", "cmsa_es_exp5.csv", "new_runs.json", "sa", "min")
# setwd("/Users/feldt/dev/BlackBoxOptim.jl/spikes/experiments")

num_params <- as.integer(args[1]);
#ranges <- matrix(eval(parse(text = args[2])), byrow=TRUE, nrow=num_params);
num_new_runs <- as.integer(args[2]);

response_column <- as.integer(args[3]); # Response that should be minimized
if(response_column < 0) {
  invert_response <- TRUE;
  response_column <- -response_column;
} else {
  invert_response <- FALSE;
}

basemax <- as.integer(args[4]); # basemax value, params larger than this are binary/categorical
csvfile <- args[5];
outfile <- args[6];

if(length(args) >= 7) {
  perform_sa <- args[7]
} else {
  perform_sa <- "not"
}

if(length(args) >= 8) {
  selection_scheme <- tolower(args[8])
} else {
  selection_scheme <- "ei"
}

# Create a result list where we will save results
result = list(analysis_date = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))


#####################################################################
# 2. Read input if the file exist
#####################################################################
csv_exists <- file.exists(csvfile)
if(csv_exists) {
  # Read input csv
  runs <- read.csv(csvfile, header = TRUE)

  num_runs <- nrow(runs)

  # Extract the X design matrix. Its columns must be the first num_params
  # columns of the runs data frame.
  X = as.matrix(runs[,1:num_params])

  # Extract the Z response vector.
  Z = as.matrix(runs[,response_column])

  # For using EI improv below we need to minimize. The sign of the response 
  # column above indicated if we need to invert.
  if(invert_response) {
    Z = -Z
  }

}


#####################################################################
# 3. Perform sensitivity analysis if requested
#####################################################################
if(csv_exists && (tolower(perform_sa) == "sa" || tolower(perform_sa) == "true")) {

  # Set general rect, mode and shape params
  #rect <- t(apply(X, 2, range, na.rm=TRUE))
  rect <- matrix(rep(0.0, 2*num_params), nrow=num_params)
  rect[,2] <- 1.0
  #mode <- apply(X , 2, mean, na.rm=TRUE)
  #mode <- rep(0.5, num_params)
  shape = rep.int(2, num_params)

  # Update rect and shape for the categorical/binary vars.
  #shape[(basemax+1):num_params] <- 0
  #rect[(basemax+1):num_params,1] <- 0
  #rect[(basemax+1):num_params,2] <- 1
  #mode[(basemax+1):num_params] <- 0.5

  cat("Performing sensitivity analysis\n")
  Ngrid <- 100
  s <- suppressWarnings(sens(X=X, Z=Z, nn.lhs=300, model=btgp, 
#    shape = shape, rect = rect,
    ngrid=Ngrid, span=0.3, BTE=c(5000,10000,10)))

  #png('main_effects_per_var.png', width = 1600, height = 1200)
  pdf('main_effects_per_var.pdf')
  plot(s, layout="sens", main="", maineff=t(1:num_params))
  dev.off()

  #png('sensitivty_per_var.png', width = 1600, height = 1200)
  pdf('sensitivity_per_var.pdf')
  plot(s, layout="sens", maineff=FALSE)
  dev.off()

  # Save sensitivity indices stats in result
  result$sa_mean_1st_order_sens_indices = colMeans(s$sens$S)
  result$sa_mean_total_sens_indices = colMeans(s$sens$T)
  result$sa_sd_1st_order_sens_indices = apply(s$sens$S, 2, sd)
  result$sa_sd_total_sens_indices = apply(s$sens$T, 2, sd)

  # Find the values for the 5% best quantiles for each param
  quantile <- 0.05
  if(invert_response == TRUE) {
    cutoffprob <- quantile;
  } else {
   cutoffprob <- 1.0 - quantile;
  }
  num_in_5 <- round(quantile*Ngrid);
  best_sa <- matrix(rep.int(0, num_in_5*num_params), nrow=num_in_5)
  for(pindex in 1:num_params) {
    q <- as.double(quantile(s$sens$ZZ.mean[,pindex], probs=c(cutoffprob)))
    # Bug below: Only saves first value instead of all of them...
    best_sa[,pindex] <- s$sens$Xgrid[which(s$sens$ZZ.mean[,pindex] <= q),pindex]
  }
  result$best_sa <- best_sa
}

#####################################################################
# 4. Build tgp model if there are any existing runs, if not we
#    generate initial points via lhs. Save points to out file.
#####################################################################
if(csv_exists && nrow(runs) > 0) {

  # Create a LHS sample of points to select from
  cat("Sampling many points to select from\n")
  Xcandidates = sample_param_space(num_params, lhs_sample_size);

  if(selection_scheme == "min") {

    XX <- Xcandidates;

  } else {

    # Build model for selecting new points
    cat("Building model for selecting points\n")
    model <- btgp(X=X, Z=Z, pred.n=FALSE, basemax = basemax, R=num_tgp_runs)
    #model2 <- btgpllm(X=X, Z=Z, pred.n=FALSE, basemax = basemax, R=num_tgp_runs)

    cat("Select a subset of points with most design value\n")
    num_points <- max(num_selected_design_points, num_new_runs);
    XX <- tgp.design(num_points, Xcandidates, model);

  }

  # Now predict in those points
  pmodel <- btgp(X=X, Z=Z, XX=XX, basemax = basemax, corr="exp", improv=TRUE,
    R=num_tgp_runs, krige = FALSE)
  #pmodel2 <- btgpllm(X=X, Z=Z, XX=XX, basemax = basemax, corr="exp", improv=TRUE,
  #  Ds2x=TRUE, R=num_tgp_runs, krige = FALSE)

  # Write the tgp tree to file
  pdf('tgp_tree.pdf')
  tgp.trees(pmodel)
  dev.off()

  # Write the posterior predictive surface per parameter to disk.
  for(i in 1:num_params) {
    name <- names(runs)[i+num_params]
    pdf(paste('posterior_', i, '_', name, '.pdf', sep=""))
    plot(pmodel, main=name, proj=c(i))
    dev.off()
  }

  # Predicted values
  ZZ = pmodel$ZZ.mean
  #ZZ2 = pmodel2$ZZ.mean

  if(selection_scheme == "min") {

    cat("Selection scheme: minimization of response\n")

    # Find the num_new_points points with minimum values and select those points.
    mins = sort(ZZ)
    index = which(ZZ %in% mins[1:num_new_runs])
    #mins2 = sort(ZZ2)
    #indmins2 = which(ZZ2 %in% mins2[1:num_new_runs])
    XXsel <- XX[index,]

  } else if(selection_scheme == "alc") {

    cat("Selection scheme: Active Learning Cohn")

  } else {

    # Default scheme is "EI", i.e. minimization according to expected improvement
    cat("Selection scheme: Expected Improvement\n")

    # Find the first ranked design point according to EI
    index <- which(pmodel$improv$rank <= num_new_runs)
    XXsel <- XX[index[1:num_new_runs], ]

  }

  # Print some info
  cat("XXsel = ", XXsel, "\n");
  cat("Predicted values for XXsel, ZZ = ", ZZ[index], "\n");

  # Print some stats for some subsets of the top list predicted
  num_min = 5
  mins = sort(ZZ)
  while(num_min <= 50) {
    index = which(ZZ %in% mins[1:num_min]);
    dvs <- XX[index,1:num_params];
    means <- colMeans(dvs);
    sds <- apply(dvs, 2, sd);
    cat("Top ", num_min, ": y = ", mean(ZZ[index]), ", ds = ", means, ", sd = ", sds, "\n");

    if(num_min == 5) {
      result$top5_response <- mean(ZZ[index]);
      result$top5_mean <- means;
      result$top5_sd <- sds;
      result$top5_min <- apply(dvs, 2, min);
      result$top5_max <- apply(dvs, 2, max);
    }

    num_min <- num_min + 5;
  }

} else {

  XXsel = sample_param_space(num_params, num_params+1);

}

result$num_rows = num_new_runs;
result$num_cols = num_params;
result$best = XXsel;

#save_matrix_to_json_file(XXsel, outfile);
save_list_to_json_file(result, outfile);