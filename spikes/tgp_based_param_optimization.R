# We will use Treed Gaussian Processes to optimize the parameter settings
# of a Compass Search implemented with GSS.
library(tgp)

#####################################################################
# 1. Optimize single run for highest success ratio with fewest fevals
#####################################################################

# Set some constants.
n <- 8 # Start with a simple 2-D cigar problem on [-10, 10]^n => diameter = 20
diameter <- 30 - (-30)
problem <- "Rosenbrock"

# 1a. First define the parameters and their ranges
num_params = 6

ranges <- matrix(c(
    # -2,  1,  # max_seconds          = 10^[-2, 1]*n    
     0,   3,  # lambda               = 4 * [2 (1), n (2), n*n (3)]
     0,   3,  # mu                   = lambda / [2 (1), 4 (2), 8 (3)]
     0,   2,  # sampler              = EigenCovarSampler (1) or CholeskyCovarSampler (2)
     0,   2,  # utilitiesFunc        = log (1) or linear (2)
     0.6, 1,  # covar_learning_rate  = [0, 1]
    -2,   0   # sigma                = 10^[-2, 0] * sample_space_diameter
), byrow=TRUE, nrow=num_params);

# 1b. Latin Hypercube sampling to get initial points
library(lhs)

sample_param_space <- function(ranges, num_samples) {
  delta <- ranges[,2] - ranges[,1]
  s <- improvedLHS(nrow(ranges), num_samples, 100)
  t(ranges[,1] + delta * s)
}

X = sample_param_space(ranges, num_params+1)

save_matrix_to_file <- function(matrix, filepath) {
  fileConn<-file(filepath)
  writeLines(toJSON(matrix), fileConn)
  close(fileConn)
}

save_matrix_to_file(X, "/Users/feldt/dev/BlackBoxOptim.jl/spikes/experiments/params.txt")

#            [,1]      [,2]      [,3]      [,4]      [,5]        [,6]
# [1,] 1.18961584 0.3892732 0.5929056 1.2097220 0.9306928 -0.02241825
# [2,] 0.05357900 0.7919653 0.9507112 1.4970144 0.8184783 -0.12791662
# [3,] 1.17934940 2.4306748 0.1344748 0.3689288 0.9779468 -0.92497880
# [4,] 0.31565622 2.1102317 1.9776501 1.1412879 0.7995758 -1.56227778
# [5,] 0.30185745 0.6737901 0.8119760 1.2950750 0.8910373 -0.30685010
# [6,] 0.01982635 2.4366438 1.1793706 0.5494892 0.9981388 -1.08508529
# [7,] 1.28083038 0.1430397 0.4893044 1.2193957 0.8905398 -0.29996725

Z = c(0.3717598530202772,8.1359371709471,7.951338981802012e-8,17405.398023972728,4431.812316119974,9.938292408687424e-8,21297.166483066914)

build_model_then_select_new_points <- function(X, Z, num_new_points = 2, reruns = 5) {

  # Now build the tgp model!
  m1 <- btgp(X=X, Z = Z, pred.n=FALSE, R=5)

  # And select the best 100 points out of a 1000 ones
  Xcandidates = sample_param_space(ranges, 1000)
  XX <- tgp.design(30, Xcandidates, m1)

  # Now predict in those points
  m1p <- btgp(X=X, Z=Z, XX=XX, corr="exp", improv=TRUE,
    Ds2x=TRUE, R=5, verb=0)

  # Predicted values
  ZZ = m1p$ZZ.mean

  # Find the 2 points with minimum values and select those points.
  mins = sort(ZZ)
  indmins = which(ZZ %in% mins[1:num_new_points])
  XXsel <- XX[indmins,]

  list(model = m1p, XXsel = XXsel)
}

# round 2:
r <- build_model_then_select_new_points(X, Z, 2, 5)
save_matrix_to_file(r$XXsel, "/Users/feldt/dev/BlackBoxOptim.jl/spikes/experiments/params2.txt")
Zsel = c(904813.0181147617,9.701535337325367e-8)
Z = c(Z, Zsel)
X = rbind(X, r$XXsel)

# round 3:
r <- build_model_then_select_new_points(X, Z, 2, 5)
save_matrix_to_file(r$XXsel, "/Users/feldt/dev/BlackBoxOptim.jl/spikes/experiments/params2.txt")
Zsel = c(0.014530563561333027,668474.6555622832)
Z = c(Z, Zsel)
X = rbind(X, r$XXsel)

# round 3:
r <- build_model_then_select_new_points(X, Z, 2, 5)
save_matrix_to_file(r$XXsel, "/Users/feldt/dev/BlackBoxOptim.jl/spikes/experiments/params2.txt")
Zsel = c(9.898542095571683e-8,9.941423547251608e-8)
Z = c(Z, Zsel)
X = rbind(X, r$XXsel)

# round 4:
r <- build_model_then_select_new_points(X, Z, 2, 5)
save_matrix_to_file(r$XXsel, "/Users/feldt/dev/BlackBoxOptim.jl/spikes/experiments/params2.txt")
Zsel = c(9.9328695150376e-8,4.662224992081074)
Z = c(Z, Zsel)
X = rbind(X, r$XXsel)

# round 5:
r <- build_model_then_select_new_points(X, Z, 2, 5)
save_matrix_to_file(r$XXsel, "/Users/feldt/dev/BlackBoxOptim.jl/spikes/experiments/params2.txt")
Zsel = c(9.9328695150376e-8,4.662224992081074)
Z = c(Z, Zsel)
X = rbind(X, r$XXsel)

# round 6:
r <- build_model_then_select_new_points(X, Z, 2, 5)
save_matrix_to_file(r$XXsel, "/Users/feldt/dev/BlackBoxOptim.jl/spikes/experiments/params2.txt")
Zsel = c(3.7956844680606583,9.98245302853361e-8)
Z = c(Z, Zsel)
X = rbind(X, r$XXsel)

# round 7:
r <- build_model_then_select_new_points(X, Z, 2, 5)
save_matrix_to_file(r$XXsel, "/Users/feldt/dev/BlackBoxOptim.jl/spikes/experiments/params2.txt")
Zsel = c(0.1794436927014203,9.219486830374518e-8)
Z = c(Z, Zsel)
X = rbind(X, r$XXsel)

# round 8:
r <- build_model_then_select_new_points(X, Z, 2, 5)
save_matrix_to_file(r$XXsel, "/Users/feldt/dev/BlackBoxOptim.jl/spikes/experiments/params2.txt")
Zsel = c(8.715664399647491e-8,8.961166930707209e-8)
Z = c(Z, Zsel)
X = rbind(X, r$XXsel)

# round 9:
r <- build_model_then_select_new_points(X, Z, 2, 5)
save_matrix_to_file(r$XXsel, "/Users/feldt/dev/BlackBoxOptim.jl/spikes/experiments/params2.txt")
Zsel = c(8.444179961161977e-8,9.967467268892997e-8)
Z = c(Z, Zsel)
X = rbind(X, r$XXsel)

# round 10:
r <- build_model_then_select_new_points(X, Z, 10, 5)
save_matrix_to_file(r$XXsel, "/Users/feldt/dev/BlackBoxOptim.jl/spikes/experiments/params2.txt")
Zsel = c(9.463207984070908e-8,2.7041464563904,5.407635734491258,318.24872395468975,118.88942110418073,9.969438777267423e-8,1.1178740598824133,0.04542230512842288,9.952328340205899e-8,9.443712033633461e-8)
Z = c(Z, Zsel)
X = rbind(X, r$XXsel)

# Do a sensitivity analysis
s <- suppressWarnings(sens(X=X, Z=Z, nn.lhs=300, model=btgp,
                           ngrid=100, span=0.3, BTE=c(5000,10000,10)))
plot(s, layout="sens", ylab="Ozone", main="main effects")
plot(s, layout="sens", ylab="Ozone", main="", maineff=t(1:6))

# round 11:
r <- build_model_then_select_new_points(X, Z, 20, 5)
save_matrix_to_file(r$XXsel, "/Users/feldt/dev/BlackBoxOptim.jl/spikes/experiments/params2.txt")
Zsel = c(757.3268417765296,166.30794883919495,120.59390128793476,45.74771951955749,197.2810853377193,795.0204446351319,4.1828766313481495e7,8.671546302426113e-8,0.681147633942023,2.4587040161001498e8,75234.7750628955,111.21589490173636,193.25840618604425,3.9858877695994885,6.80307951814857e-8,0.885256898154827,6.08163027208339e-8,3.9858877695994885,855.0580725567552,8.866038734145751e-8)
Z = c(Z, Zsel)
X = rbind(X, r$XXsel)

# round 12:
r <- build_model_then_select_new_points(X, Z, 20, 5)
save_matrix_to_file(r$XXsel, "/Users/feldt/dev/BlackBoxOptim.jl/spikes/experiments/params2.txt")
Zsel = c(13.526925584549623,6.066077506163768e-8,576118.564677597,0.002177707848724537,293.45369510536506,7.226549532422246e-8,3.9858877695994885,9.919971399533284e-8,8.286270534063318e-8,9.591188463849613e-8,54.85021865375828,7928.029585030435,3.985977991369322,3.9858877695994885,247658.93534100323,246632.06705601697,9.710133673960906e-8,49.448993999883086,5.138210345705269,6.963342285746061e-5)
Z = c(Z, Zsel)
X = rbind(X, r$XXsel)

# There seems to be a lot of variation though since the sensitivity models
# vary a lot between different runs.
# Lets run a larger sensitivity analysis:
s <- suppressWarnings(sens(X=X, Z=Z, nn.lhs=600, model=btgp,
                           ngrid=200, span=0.3, BTE=c(5000,10000,10)))
plot(s, layout="sens", ylab="Ozone", main="main effects")
plot(s, layout="sens", ylab="Ozone", main="", maineff=t(1:6))

# Conclusions:
#   1. covar_learning_rate needs to be at least 0.90, but lower responses if high
#   2. sigma should be at least -1, but lower responses all the way up to 0
#   3. utility func should be log
#   4. sampler should be cholesky
#   5. lambda should be higher rather than lower
#   6. the mu divisor should be higher rather than lower but least sensitive to this one
#
# So lets fix sampler and utility func according to the above conclusions and 
# do a new run with parameters from
num_params = 6

ranges <- matrix(c(
    # -2,  1,  # max_seconds          = 10^[-2, 1]*n    
     0,   3,  # lambda               = 4 * [2 (1), n (2), n*n (3)]
     0,   3,  # mu                   = lambda / [2 (1), 4 (2), 8 (3)]
     0,   2,  # sampler              = EigenCovarSampler (1) or CholeskyCovarSampler (2)
     0,   2,  # utilitiesFunc        = log (1) or linear (2)
     0.6, 1,  # covar_learning_rate  = [0, 1]
    -2,   0   # sigma                = 10^[-2, 0] * sample_space_diameter
), byrow=TRUE, nrow=num_params);


# GP paket att testa: tgp, dynatree, GPfit, mlegp