## Monte Carlo simulation for estimating the distribution of the P statistic
## (16), necessary for computing the Quantile Universal Threshold.
## (10'000 realizations distributed in 10 job arrays on the cluster,
## use aggregate.qutMC.R to aggregate results)


# SLURM array number:
k <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")) # k between 1 and 10

# parameters:
type <- Sys.getenv("type") # in batch script: --export = type="smallGauss", ALL e.g.
MCrep <- 1000

# necessary packages and functions
library(lpSolve)
library(doParallel)
registerDoParallel(cores = detectCores())
source("generateX.R")
BP <- function(X, y) {
    p <- ncol(X)
    out.lp <- lpSolve::lp(objective.in = rep(1, 2*p), const.mat = cbind(X, -X), 
                          const.dir = "==", const.rhs =y )
    if (out.lp$status == 2) stop("linear system y = X beta + G gamma has no solution; choose larger q (e.g. q >= nrow(X))" )
    sol <- out.lp$solution[1:p] - out.lp$solution[-(1:p)]
    sol
}

# design matrix:
X <- generateX(type = type)
n <- nrow(X)
p <- ncol(X)

# lass0 parameters:
intercept <-  TRUE # ensures that the noise dictionaries are also mean-centered - even though X already is
standardizeX <- FALSE # already done by generateX when necessary
q <- n
M <- 30
if (type == "TV300") {
    source("comp.TV300.Gnorm.R")
    N <- 1.e4
    standardizeG <- comp.TV300.Gnorm(N = N, q = q, alpha = 0.05, 
                                     parallel = TRUE, seeds = 1:N) # returns the desired Euclidean norm for columns of G
} else {
    standardizeG <- TRUE
}

# Monte Carlo simulation:
# (equiv. to calling:
# MCoutput <- qut.MC(X, q = q, intercept = intercept, 
#                   standardizeX = standardizeX, standardizeG = standardizeG, 
#                   MCrep = MCrep, GEVapprox = FALSE, parallel = TRUE)
# with qut.MC function from lass0 package)
allMC <- foreach(r = 1:MCrep, .combine = "c", .packages = "lpSolve") %dopar% {
    print(paste0("MC rep", r))
    set.seed((k-1) * MCrep + r)
    eps <- rnorm(n)
    if (intercept) eps <- eps - mean(eps)
    BPsols <- array(NA, dim=c(p+q, M))
    for (m in 1:M) {
        G <- matrix(rnorm(n*q), n, q)
        if (intercept) G <- t(t(G) - colMeans(G))
        if (is.numeric(standardizeG)) {
            G <- standardizeG * t(t(G) / apply(G, 2, function(v) sqrt(sum(v^2))))
        } else if (standardizeG) {
            G <- t(t(G) / apply(G, 2, sd))
        }
        BPsols[, m] <- BP(cbind(X, G), eps)
    }
    Betas <- BPsols[1:p, ]
    medBeta <- apply(Betas, 1, median)
    Gammas <- BPsols[-(1:p), ]
    gammas <- as.vector(Gammas)
    madGammas <- mad(gammas[gammas != 0])
    max(abs(medBeta))/madGammas
}

MCoutput <- list()
MCoutput$allMC <- allMC
MCoutput$lass0setting <- list(q = q, M = M, sigma = NULL, intercept = intercept, 
                              standardizeX = standardizeX, standardizeG = standardizeG)



# save results
assign(paste0("qutMC.", type, ".array", k), MCoutput)
save(list=paste0("qutMC.", type, ".array", k), 
     file=paste0("qutMC.", type, ".array", k, ".Rda"))
