library(optparse, quietly = TRUE, warn.conflicts = FALSE, verbose = FALSE)
library(plyr, quietly = TRUE, warn.conflicts = FALSE, verbose = FALSE)

option_list <- list(

  make_option(c("-e", "--experimentname"), type="character", default="exp",
    help="Name of experiment to be analysed, possibly with path"),

  make_option(c("-s", "--inputseparator"), type="character", default=",",
    help="Separator to use when reading in data"),

  make_option(c("--details"), action="store_true",
    help="Give detailed info per problem"),

  make_option(c("-t", "--ftol"), type="double", default=1e-7,
    help="Fitness tolerance used as target in optimizations")

)

args <- parse_args(OptionParser(option_list = option_list))

give_details <- !is.null(args$details) 

# Create the two csv file names
summaryfile = paste(args$experimentname, "_summary.csv", sep="")
histfile = paste(args$experimentname, "_runs.csv", sep="")


# Read summary data.
cat("Reading summary input csv", args$inputfile);
time_read <- system.time( dsummary <- read.csv(summaryfile, header = TRUE, 
  sep=args$inputseparator) );
cat(" (", time_read[3], " seconds)\n", sep = "");

num_runs <- nrow(dsummary)
cat("Number of runs in experiment: ", num_runs, "\n")

problems <- unique(dsummary$Problem)
num_problems <- length(problems)
cat("Number of problems in experiment: ", num_problems, "\n")
cat("Problems: ", paste(problems), "\n")

dims <- unique(dsummary$Dimension)
num_dims <- length(dims)
cat("Number of dimensions in experiment: ", num_dims, "\n")
cat("Dimensions: ", paste(dims), "\n")


# Read fitness history data.
cat("\nReading fitness history input csv", args$inputfile);
time_read <- system.time( dhist <- read.csv(histfile, header = TRUE, 
  sep=args$inputseparator) );
cat(" (", time_read[3], " seconds)\n", sep = "");
cat("Number of entries in history =", nrow(dhist), "\n")


#####################################################################
# Aggregate median, mean, std, min and max fitness as well as
# success rate and execution time per problem and dim.
#####################################################################

ratio_fitness_below_ftol <- function(fitnesses, ftol = args$ftol) {
  sum(fitnesses < ftol) / length(fitnesses)
}

detailed_summary_by_problem_and_dim <- ddply(dsummary, c("Problem", "Dimension"), summarise,
               SuccessRate   = 100.0*ratio_fitness_below_ftol(Fitness),
               MedianFitness = median(Fitness),
               MeanFitness   = mean(Fitness),
               StdDevFitness = sd(Fitness),
               MinFitness    = min(Fitness),
               MaxFitness    = max(Fitness),
               MedianFevals  = median(FuncEvals),
               MeanFevals    = mean(FuncEvals),
               StdDevFevals  = sd(FuncEvals),
               MinFevals     = min(FuncEvals),
               MaxFevals     = max(FuncEvals),
               MedianElapsedTime = median(ElapsedTime),
               NumReps       = length(Fitness))

if(give_details) {
  print.data.frame(detailed_summary_by_problem_and_dim[,c("Problem", "Dimension", 
    "SuccessRate", "MedianFitness", "MedianFevals", "MedianElapsedTime", "NumReps")]);
}

detailed_summary_by_dim <- ddply(dsummary, c("Dimension"), summarise,
               SuccessRate   = 100.0*ratio_fitness_below_ftol(Fitness),
               MedianFitness = median(Fitness),
               MeanFitness   = mean(Fitness),
               StdDevFitness = sd(Fitness),
               MinFitness    = min(Fitness),
               MaxFitness    = max(Fitness),
               MedianFevals  = median(FuncEvals),
               MeanFevals    = mean(FuncEvals),
               StdDevFevals  = sd(FuncEvals),
               MinFevals     = min(FuncEvals),
               MaxFevals     = max(FuncEvals),
               MedianElapsedTime = median(ElapsedTime),
               NumReps       = length(Fitness))

print.data.frame(detailed_summary_by_dim[,c("Dimension", "SuccessRate", "MedianFitness", "MedianFevals", "MedianElapsedTime", "NumReps")])
