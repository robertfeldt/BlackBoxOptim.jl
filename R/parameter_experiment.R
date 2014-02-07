# Invoke tgp to do Treed Gaussian Process regression for prediction and design
# of parameter experiments on black box optimizers.
# Copyright (c) 2013-2014 Robert Feldt, robert.feldt@gmail.com
library(tgp)
library(rjson)
library(lhs)
library(optparse, quietly = TRUE, warn.conflicts = FALSE, verbose = FALSE)
library(fields)

#####################################################################
# 0. Define constants and functions
#####################################################################
num_tgp_runs <- 3                 # Number of runs when building tgp models
improvedlhs_num_repeats <- 5      # Number of points to consider for each added design point with improvedLHS
lhs_sample_size <- 500            # Number of candidate points sampled when selecting design points
oversampling_factor_for_minquick <- 10 # Number of times more candidate points to include when minimizing, from which the quickest is then selected
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
# args <- c("6", "1", "-17", "6", "cmsa_es_exp9.csv", "new_runs.json", "sa", "min")
# args <- c("8", "1", "17", "8", "cmsa_es_Griewank_32_8_params.csv", "new_runs.json", "sa", "min")
# setwd("/Users/feldt/feldt/research/projects/optimizing_hyperparameters_via_treed_gaussian_processes/experiments/1_5problems_6_dimensions")

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

if(length(args) >= 9) {
  response_time_column <- tolower(args[9])
} else {
  response_time_column <- 2*num_params+5-1
}

if(selection_scheme == "ei") {
  improv_flag = TRUE;
  alc_flag = FALSE;
} else if(selection_scheme == "alc") {
  alc_flag = TRUE;
  improv_flag = FALSE;
} else {
  alc_flag = FALSE;
  improv_flag = FALSE;  
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
  colnames(X) <- names(runs)[(num_params+1):(2*num_params)]

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
  s <- suppressWarnings(sens(X=X, Z=Z, nn.lhs=(100*num_params), model=btgp, 
#    shape = shape, rect = rect,
    ngrid=Ngrid, span=0.3, BTE=c(5000,10000,10)))

  # Save sensitivity indices stats in result
  result$sa_mean_1st_order_sens_indices = colMeans(s$sens$S)
  result$sa_mean_total_sens_indices = colMeans(s$sens$T)
  result$sa_sd_1st_order_sens_indices = apply(s$sens$S, 2, sd)
  result$sa_sd_total_sens_indices = apply(s$sens$T, 2, sd)
  result$sa_var_order_1st <- (1:num_params)[order(result$sa_mean_1st_order_sens_indices, decreasing = TRUE)]
  result$sa_var_order_total <- (1:num_params)[order(result$sa_mean_total_sens_indices, decreasing = TRUE)]

  # Save the indices to the 4 most sensitive params (according to main effect)
  most_sensitive <- result$sa_var_order_1st[1:4]

  pdf('main_effects_per_var.pdf')
  plot(s, layout="sens", main="", maineff=t(most_sensitive))
  dev.off()

  pdf('sensitivity_per_var.pdf')
  plot(s, layout="sens", maineff=FALSE)
  dev.off()

  # Find the values for the best quantiles for each param and save in results.
  # Since we have made sure to invert Z values above if we are maximizing we are
  # always minimizing here and thus should select the bottom quantile.
  quantile <- 0.05
  num_in_top <- round(quantile*Ngrid);
  best_sa <- matrix(rep.int(0, num_in_top*num_params), nrow=num_in_top)
  for(pindex in 1:num_params) {
    q <- as.double(quantile(s$sens$ZZ.mean[,pindex], probs=c(quantile)))
    indices_that_give_best_outcome <- which(s$sens$ZZ.mean[,pindex] <= q)
    perm <- order(s$sens$ZZ.mean[indices_that_give_best_outcome,pindex])
    top_indices <- indices_that_give_best_outcome[perm][1:num_in_top]
    best_sa[,pindex] <- s$sens$Xgrid[top_indices,pindex]
  }
  result$best_sa <- best_sa
  result$best_sa_num_rows <- num_in_top
  result$best_sa_num_cols <- num_params
}

#####################################################################
# 4. Build tgp model if there are any existing runs, if not we
#    generate initial points via lhs. Save points to out file.
#####################################################################
if(csv_exists && nrow(runs) > 0) {

  # Create a LHS sample of points to select from
  cat("Sampling many points to select from\n")
  Xcandidates = sample_param_space(num_params, lhs_sample_size, improvedlhs_num_repeats);

  if(selection_scheme == "min" || selection_scheme == "minquick") {

    XX <- Xcandidates;

  } else {

    # Build model for selecting new points
    cat("Building model for selecting points\n")
    model <- btgp(X=X, Z=Z, pred.n=FALSE, basemax = basemax, R=num_tgp_runs)

    cat("Select a subset of points with most design value\n")
    num_points <- max(num_selected_design_points, num_new_runs);
    XX <- tgp.design(num_points, Xcandidates, model);

  }

  # Now predict in those points
  pmodel <- btgp(X=X, Z=Z, XX=XX, basemax = basemax, corr="exp", 
    improv = improv_flag, Ds2x = alc_flag, R=num_tgp_runs, krige = FALSE)

  # Write the tgp tree to file
  #pdf('tgp_tree.pdf')
  #tgp.trees(pmodel)
  #dev.off()

  # Write the posterior predictive surface for the main effect to 2nd and from
  # main effect to 3rd, if we performed a sensitivity analysis above.
  if("sa_var_order_1st" %in% names(result)) {
    v1 <- result$sa_var_order_1st[1]
    namev1 <- names(runs)[v1+num_params]
    best_values <- result$best_sa[1,]

    for(i in 2:num_params) {
      vi <- result$sa_var_order_1st[i]
      namevi <- names(runs)[vi+num_params]

      # Plot posterior surface plot
      pdf(paste('posterior_surface_1_', i, '.pdf', sep=""))
      plot(pmodel, main=paste(substr(namev1,1,6), " vs. ", substr(namevi, 1,6), sep=""), proj=c(v1, vi))
      dev.off()

      # Plot posterior image 2d plot when fixing the "other" vars at their best values
      fixed_vars <- setdiff(1:num_params, c(v1, vi))
      x1 = x2 = seq(0.0, 1.0, length.out=100)
      grid <- expand.grid(x = x1, y = x2)
      gridxx = matrix(rep(0.0, num_params*nrow(grid)), nrow=nrow(grid))
      gridxx[,v1] <- grid$x
      gridxx[,vi] <- grid$y
      for(j in 1:length(fixed_vars)) {
        gridxx[,fixed_vars[j]] <- rep(best_values[fixed_vars[j]], nrow(grid))
      }
      o <- predict(pmodel, XX=gridxx, pred.n=FALSE)
      filename = paste('posterior_image_1_', i, '.pdf', sep="")
      pdf(filename)
      image.plot(x1, x2, matrix(o$ZZ.mean, nrow=length(x1)))
      dev.off()
    }
  }

  # Predicted values
  ZZ = pmodel$ZZ.mean

  if(selection_scheme == "min") {

    cat("Selection scheme: minimization of response\n")

    # Find the num_new_points points with minimum values and select those points.
    mins = sort(ZZ)
    index = which(ZZ %in% mins[1:num_new_runs])
    XXsel <- XX[index,]

  } else if(selection_scheme == "minquick") {

    cat("Selection scheme: quickest candidate when minimizing response\n")

    # Find the num_new_points*oversampling_factor_for_minquick points with 
    # minimum values and select those points.
    mins = sort(ZZ)
    index_min = which(ZZ %in% mins[1:(num_new_runs*oversampling_factor_for_minquick)])

    # Extract the response time column
    Ztime = as.matrix(runs[,response_time_column])

    # Now build a model for the median response time and then predict the response
    # time for the selected candidate set.
    cat("Building model of response time\n")
    tmodel <- btgp(X=X, Z=Ztime, XX=XX, basemax = basemax, corr="exp", 
      improv = FALSE, Ds2x = FALSE, R=num_tgp_runs, krige = FALSE)

    # Get the predicted times
    pred.times = tmodel$ZZ.mean

    cat("Oversampled set:\n"); print.table(XX[index_min,]);
    cat("Predicted response times for oversampled set: ", pred.times[index_min], "\n")

    # Find the point(s) among them with minimum predicted time
    min_times = sort(pred.times[index_min])
    index = which(pred.times %in% min_times[1:num_new_runs])
    result$predicted_times <- min_times[1:num_new_runs]
    XXsel <- XX[index,]

  } else if(selection_scheme == "alc") {

    cat("Selection scheme: Active Learning Cohn")

    alc_sorted <- sort(pmodel$Ds2x, decreasing = FALSE)
    min_to_include <- alc_sorted[num_new_runs]
    index <- which(pmodel$Ds2x <= min_to_include)
    XXsel <- XX[index, ]

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
  result$predicted_responses <- ZZ[index]
  if(selection_scheme == "minquick") {
    cat("Predicted time XXsel, pred.times = ", pred.times[index], "\n");
  }

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

  XXsel = sample_param_space(num_params, num_new_runs, improvedlhs_num_repeats);

}

result$num_rows = num_new_runs;
result$num_cols = num_params;
result$best = XXsel;

#save_matrix_to_json_file(XXsel, outfile);
save_list_to_json_file(result, outfile);