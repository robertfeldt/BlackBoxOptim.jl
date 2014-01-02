library(optparse, quietly = TRUE, warn.conflicts = FALSE, verbose = FALSE)

option_list <- list(
  make_option(c("-e", "--experimentname"), type="character", default="exp",
    help="Name of experiment to be analysed, possibly with path"),

  make_option(c("-s", "--inputseparator"), type="character", default=",",
    help="Separator to use when reading in data")

)

args <- parse_args(OptionParser(option_list = option_list))

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

dims <- unique(dsummary$Dimensions)
num_dims <- length(dims)
cat("Number of dimensions in experiment: ", num_dims, "\n")
cat("Dimensions: ", paste(dims), "\n")


# Read fitness history data.
cat("\nReading fitness history input csv", args$inputfile);
time_read <- system.time( dhist <- read.csv(histfile, header = TRUE, 
  sep=args$inputseparator) );
cat(" (", time_read[3], " seconds)\n", sep = "");
cat("Number of entries in history =", nrow(dhist), "\n")
