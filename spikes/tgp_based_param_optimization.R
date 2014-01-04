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

# X:
#            [,1]       [,2]       [,3]        [,4]      [,5]         [,6]
 #[1,] 1.18961584 0.38927316 0.59290557 1.209722046 0.9306928 -0.022418249
 #[2,] 0.05357900 0.79196530 0.95071116 1.497014398 0.8184783 -0.127916624
 #[3,] 1.17934940 2.43067476 0.13447478 0.368928760 0.9779468 -0.924978799
 #[4,] 0.31565622 2.11023174 1.97765010 1.141287866 0.7995758 -1.562277782
 #[5,] 0.30185745 0.67379012 0.81197604 1.295075046 0.8910373 -0.306850100
 #[6,] 0.01982635 2.43664377 1.17937057 0.549489173 0.9981388 -1.085085294
 #[7,] 1.28083038 0.14303969 0.48930442 1.219395728 0.8905398 -0.299967247
 #[8,] 0.96103657 2.13781702 1.18298674 1.871317507 0.7742555 -1.994522352
 #[9,] 0.29074505 1.86268614 0.57113889 1.614156265 0.9508397 -1.017206343
#[10,] 0.30165002 1.18227717 1.80627057 0.414625693 0.8399694 -0.605599993
#[11,] 0.25635239 2.43976712 0.90751827 0.656244927 0.8042611 -0.070436746
#[12,] 0.14897294 0.60684422 1.64666838 1.071086098 0.9889911 -1.034922166
#[13,] 1.96258863 0.05557772 1.41107959 0.792606998 0.9731599 -1.531104588
#[14,] 1.01887927 0.03355568 1.13452865 0.435104786 0.9967259 -0.493864754
#[15,] 1.00766400 0.05734124 1.10698837 1.777609473 0.9165139 -1.351103326
#[16,] 1.12142113 0.81683586 1.87009143 1.227747331 0.8803141 -1.941064648
#[17,] 1.09711072 2.41143511 0.36522081 1.265987252 0.9991515 -1.944089373
#[18,] 2.46561448 0.98610919 0.15253226 1.141791265 0.9967521 -1.109170138
#[19,] 2.67186412 1.27546991 1.29160707 0.018927502 0.9282677 -1.537844805
#[20,] 1.77888766 1.13497963 0.16850829 0.441224194 0.9325872 -0.285776854
#[21,] 1.16049691 1.64065386 1.46404822 0.360929175 0.9544241 -1.760813987
#[22,] 2.96672428 0.12688306 0.88577929 0.562006249 0.8969421 -0.840486718
#[23,] 0.08407914 0.50106608 1.98358634 1.434563705 0.8045431 -1.204911770
#[24,] 2.66146702 1.29943893 1.54004337 0.041588930 0.8535306 -1.571047902
#[25,] 1.11502680 0.04665895 1.17465536 1.893624761 0.7128486 -0.656199935
#[26,] 2.37878987 0.46483416 1.89675341 1.206134287 0.7530095 -1.584467990
#[27,] 0.05683583 1.47770759 1.52677282 1.966989868 0.6852384 -0.774613261
#[28,] 1.18230036 0.06361809 1.10380480 1.778395372 0.9211261 -1.560600036
#[29,] 2.12462033 1.42561109 0.28229560 1.121890867 0.9439226 -1.615493874
#[30,] 0.41392018 2.90347567 0.99618014 1.254802855 0.9191104 -1.618116701
#[31,] 2.20032156 1.22998245 0.59117522 0.009913836 0.9971820 -0.840063805
#[32,] 0.61165358 1.44368244 0.04121816 1.181784388 0.9994378 -0.442418713
#[33,] 0.15798762 0.77164932 1.41815518 0.766936580 0.9345916 -0.818438226
#[34,] 2.36044870 0.50820745 0.08722591 1.897187500 0.7552901 -0.975065856
#[35,] 1.42422691 0.77147587 1.96769099 1.344739941 0.6271676 -0.820008787
#[36,] 1.55814072 0.02463744 1.39794097 0.607069714 0.7939946 -0.290070696
#[37,] 1.28343390 1.84057143 1.50246331 0.089554383 0.6715661 -0.252795431
#[38,] 2.61420632 1.00801136 1.39253202 0.383694040 0.6401922 -0.984781873
#[39,] 2.38952024 0.12059274 0.74003969 1.013102578 0.6768125 -0.218760597
#[40,] 2.18993152 0.88079479 0.95541085 1.930665030 0.8609671 -1.703873461
#[41,] 2.25450158 1.19602534 0.35397828 0.179818753 0.8285742 -0.179851748
#[42,] 0.32612142 1.93662356 0.61409386 1.870271392 0.8877730 -1.238452211
#[43,] 2.47021475 0.67975542 1.02367101 1.755884901 0.6142110 -1.241807858
#[44,] 2.34559803 1.98880770 1.82678945 0.983845599 0.6657926 -1.589479285
#[45,] 2.44601452 1.82160633 0.52660928 0.885459092 0.6116558 -0.065694111
#[46,] 1.10249493 0.18008418 0.61027917 1.377439734 0.9735046 -0.877331879
#[47,] 0.80542407 1.99520140 1.84147873 0.017640092 0.9160769 -1.089059685
#[48,] 1.75160008 2.96622285 0.13828634 0.696578999 0.9129867 -1.622239171
#[49,] 1.76177763 0.95874546 0.18716017 1.969700031 0.9066582 -1.174397573
#[50,] 1.00993318 2.63428224 0.52219604 0.012138859 0.9179755 -0.682994223
#[51,] 1.61632671 0.07650151 1.77472477 0.930139154 0.8978267 -1.553625460
#[52,] 0.03664195 1.91289103 0.55261161 0.806580398 0.8923274 -0.009212121
#[53,] 2.72590750 1.48192988 0.18170511 0.400970847 0.8931116 -0.885818237

# Z:
# [1] 3.717599e-01 8.135937e+00 7.951339e-08 1.740540e+04 4.431812e+03
# [6] 9.938292e-08 2.129717e+04 9.048130e+05 9.701535e-08 1.453056e-02
#[11] 6.684747e+05 9.898542e-08 9.941424e-08 9.932870e-08 4.662225e+00
#[16] 3.795684e+00 9.982453e-08 1.794437e-01 9.219487e-08 8.715664e-08
#[21] 8.961167e-08 9.989251e-08 2.762606e+05 9.463208e-08 2.704146e+00
#[26] 5.407636e+00 3.182487e+02 1.188894e+02 9.969439e-08 1.117874e+00
#[31] 4.542231e-02 9.952328e-08 9.443712e-08 7.573268e+02 1.663079e+02
#[36] 1.205939e+02 4.574772e+01 1.972811e+02 7.950204e+02 4.182877e+07
#[41] 8.671546e-08 6.811476e-01 2.458704e+08 7.523478e+04 1.112159e+02
#[46] 1.932584e+02 3.985888e+00 6.803080e-08 8.852569e-01 6.081630e-08
#[51] 3.985888e+00 8.550581e+02 8.866039e-08 1.352693e+01 6.066078e-08
#[56] 5.761186e+05 2.177708e-03 2.934537e+02 7.226550e-08 3.985888e+00
#[61] 9.919971e-08 8.286271e-08 9.591188e-08 5.485022e+01 7.928030e+03
#[66] 3.985978e+00 3.985888e+00 2.476589e+05 2.466321e+05 9.710134e-08
#[71] 4.944899e+01 5.138210e+00 6.963342e-05

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
ranges2 <- matrix(c(
     8,   (4*n*n),  # lambda         = 4 * [2 (1), n (2), n*n (3)]
     2,   8,  # mu                   = lambda / [2, 8]
     2,   2,  # sampler              = EigenCovarSampler (1) or CholeskyCovarSampler (2)
     0,   0,  # utilitiesFunc        = log (1) or linear (2)
     0.80, 1,  # covar_learning_rate 
    0.2,   3   # sigma                = [0.2, 3] * sample_space_diameter
), byrow=TRUE, nrow=num_params);

PS = sample_param_space(ranges2, 4+1)
X = PS[,c(1,2,5,6)]

# GP paket att testa: tgp, dynatree, GPfit, mlegp