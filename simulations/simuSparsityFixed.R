#### simu for fixed sparsity

rm(list=ls())

################ simulation parameters:

## sparsity
s <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
## simu type
type <- Sys.getenv("type") # "smallGauss", "wideGauss", "TV300", "riboflavin"
## design:
source("generateX.R")
set.seed(2019)
X <- generateX(type = type)
n <- nrow(X)
p <- ncol(X) 
## noise level
sigma <- 1
## number of rep
R <- 500
## targeted FDR:
alpha <- 0.05
## lass0 parameters:
M <- 30
q <- n
GEVapprox <- FALSE
ols <- FALSE
## estimators:
if (type %in% c("smallGauss", "wideGauss")) {
    amp <- 3/4
    est.names <- c( "lassoZero", 
                    "lasso",
                    "MXknockoffs.off0",
                    "MXknockoffs.off1",
                    "SLOPE.gauss",
                    "stabilitySelection",
                    "SCAD")
    # NOTE: 
    # Choice of tuning parameters:
    # lassoZero: QUT (based on pivotized statistic P (16))
    # lasso: GIC
    # stabilitySelection: 
    # adaptiveLasso: GIC
    # SCAD: GIC
    # 
    # (when GIC is used, sigma is estimated with qut package)
    # 
    # (For sign recovery: correct signs are given to true positives detected by stability selection and knockoff filter)
    # 
} else if (type == "TV300") {
    amp <- 3
    est.names <- c( "lassoZero", # QUT with pivotized statistics
                    "lasso", # note: when GIC is used, sigma is estimated with qut package
                    "SCAD", 
                    "wildBinarySeg",
                    "knockoff.Xfixed.off0")
} else if (type == "riboflavin") {
    amp <- 2
    est.names <- c("lassoZero",
                   "lasso", 
                   "SCAD",
                   "stabilitySelection")
}


#####################################




## load packages & functions
library(doParallel)
library(glmnet)
library(abind)
library(knockoff)
library(SLOPE)
library(stabs)
library(wbs)
library(ncvreg) # for SCAD and MCP
library(lass0)
library(qut)
source("lasso.byGIC.R")
source("ncv.byGIC.R")
if (type == "TV300") source("comp.TV300.Gnorm.R")

## parallel
registerDoParallel(cores=detectCores())

## support S^0 for TV300
if (type == "TV300") {
    if(s==0){
        S <- integer(0)
    }else{
        S <- round(seq(0, n, length=s + 2))[-c(1, s+2)]
    }
}


## lass0 options
standardizeX <- FALSE # already done when necessary
if (type == "TV300") {
    standardizeG <- comp.TV300.Gnorm(N = 1.e4, q = q, alpha = alpha, 
                                     parallel = TRUE, seeds = 1:1.e4)
} else{
    standardizeG <- TRUE 
}


## load qutMC realizations:
load(paste0("qutMC.", type, ".Rda"))
qut.MC.output <- get(paste0("qutMC.", type))
rm(list = paste0("qutMC.", type))

## knockoffs par
if (type == "TV300") {
    stat.used <- stat.lasso_lambdasmax
} else {
    stat.used <- stat.glmnet_coefdiff
}


## criteria
comp.crit <- function(signs.hat, signs){
    S <- which(signs != 0)
    Shat <- which(signs.hat != 0)
    if (length(S) == 0) {
        tpp <- NA
    } else {
        tpp <- sum(S %in% Shat) / length(S) 
    }
    fdp <- sum(!(Shat %in% S)) / max(1, length(Shat))
    sign.recovery <- all(signs.hat == signs)
    supp.recovery <- all((signs.hat == 0) == (signs == 0))
    modelsize <- length(Shat)
    Nfalse <- sum(!(Shat %in% S))
    return(c(sign.recovery, supp.recovery, tpp, fdp, modelsize, Nfalse))
}


## combine function
comb <- function(...){
    abind(..., along=3)
}

## simu
results <- foreach(r = 1:R, .combine="comb", .multicombine=TRUE, .packages=c("glmnet", "lpSolve", "lass0")) %dopar% {
    set.seed(r)
    print(r)
    if(type != "TV300") S <- sample(1:p, s) 
    betaS <- amp * sign(rnorm(s))
    beta <- rep(0, p)
    beta[S] <- betaS
    signs <- sign(beta)
    mu <- X[, S, drop=FALSE] %*% betaS
    y <- mu + rnorm(n, 0, sigma)
    
    if ("lassoZero" %in% est.names) {
        print(c(r, "start lasso-zero"))
        lassoZero <- lass0(X, y, alpha = alpha, q = q, M = M,
                           standardizeX = standardizeX, standardizeG = standardizeG,
                           qut.MC.output = qut.MC.output, GEVapprox = GEVapprox,
                           ols = ols)$coefficients
        lassoZero <- sign(lassoZero)
        print(c(r, "end lasso-zero"))
    }
    if (any(c("lasso", "SLOPE.gauss", "adaptiveLasso", "SCAD") %in% est.names)) { 
        # when estimation of sigma required
        print(c(r, "start sigma estimation"))
        if (type == "TV300") {
            sigma.est <- mad(diff(y))
        } else {
            sigma.est <- sigmarcv(y, X, cv = TRUE)$sigmahat # cv = TRUE -> estimated as suggested by Reid et al. 2013. 
        }
        print(c(r, "end sigma estimation"))
    }
    if("lasso" %in% est.names) {
        print(c(r, "start lasso"))
        lasso <- lasso.byGIC(X, y, sigma = sigma.est, intercept = TRUE, 
                             standardize = standardizeX)$betahat
        lasso <- sign(lasso)
        print(c(r, "end lasso"))
    }
    if ("MXknockoffs.off0" %in% est.names){
        print(c(r, "start MXknockoffs.off0"))
        off0.selected <- knockoff.filter(X, y, statistic = stat.used,
                                         fdr = alpha, offset = 0)$selected
        off0.td <- off0.selected[which(beta[off0.selected] != 0)] # true discoveries
        off0.fd <- off0.selected[which(beta[off0.selected] == 0)] # false discoveries
        MXknockoffs.off0 <- rep(0, p)
        MXknockoffs.off0[off0.td] <- sign(beta[off0.td])
        MXknockoffs.off0[off0.fd] <- 1 # (anything different from 0)
        print(c(r, "end MXknockoffs.off0"))
    }
    if ("MXknockoffs.off1" %in% est.names){
        print(c(r, "start MXknockoffs.off1"))
        off1.selected <- knockoff.filter(X, y, statistic = stat.used,
                                         fdr = alpha, offset = 1)$selected
        off1.td <- off1.selected[which(beta[off1.selected] != 0)] # true discoveries
        off1.fd <- off1.selected[which(beta[off1.selected] == 0)] # false discoveries
        MXknockoffs.off1 <- rep(0, p)
        MXknockoffs.off1[off1.td] <- sign(beta[off1.td])
        MXknockoffs.off1[off1.fd] <- 1 
        print(c(r, "end MXknockoffs.off1"))
    }
    if ("SLOPE.gauss" %in% est.names) {
        print(c(r, "start SLOPE.gauss"))
        SLOPE.gauss <- SLOPE(X, y, fdr = alpha, lambda = "gaussian", 
                             sigma = sigma.est, normalize = TRUE)$beta
        SLOPE.gauss <- sign(SLOPE.gauss)
        print(c(r, "end SLOPE.gauss"))
    }
    if("stabilitySelection" %in% est.names){
        print(c(r, "start stabilitySelection"))
        stab.selected <- stabsel(as.matrix(X), y, cutoff = 0.6, PFER = 1)$selected
        stab.td <- stab.selected[which(beta[stab.selected] != 0)]
        stab.fd <- stab.selected[which(beta[stab.selected] == 0)]
        stabilitySelection <- rep(0, p)
        stabilitySelection[stab.td] <- sign(beta[stab.td])
        stabilitySelection[stab.fd] <- 1 
        print(c(r, "end stabilitySelection"))
    }
    if ("SCAD" %in% est.names) {
        SCAD <- sign(ncv.byGIC(X, y, penalty = "SCAD", sigma = sigma.est, 
                               standardize = standardizeX))
    }
    if ("wildBinarySeg" %in% est.names){
        print(c(r, "start wild binary seg"))
        wild.res <- wbs(y)
        wbs.cp <- changepoints(wild.res)
        wbs.cp <- sort(wbs.cp$cpt.th[[1]])
        wbs.td <- wbs.cp[which(beta[wbs.cp] != 0)]
        wbs.fd <- wbs.cp[which(beta[wbs.cp] == 0)]
        wildBinarySeg <- rep(0, p)
        wildBinarySeg[wbs.td] <- sign(beta[wbs.td])
        wildBinarySeg[wbs.fd] <- 1 
        print(c(r, "end wild binary seg"))
    }
    if("knockoff.Xfixed.off0" %in% est.names){
        knockoff.Xfixed.off0 <- knockoff.filter(X, y, knockoffs = create.fixed,
                                           statistic = stat.used,
                                           fdr = alpha, offset = 0)$selected
        
        
        Xfixed.off0.td <- knockoff.Xfixed.off0[which(beta[knockoff.Xfixed.off0] != 0)] # true discoveries
        Xfixed.off0.fd <- knockoff.Xfixed.off0[which(beta[knockoff.Xfixed.off0] == 0)] # false discoveries
        knockoff.Xfixed.off0 <- rep(0, p)
        knockoff.Xfixed.off0[Xfixed.off0.td] <- sign(beta[Xfixed.off0.td])
        knockoff.Xfixed.off0[Xfixed.off0.fd] <- 1 # (anything different from 0)
    }
    
    all.est <- do.call(list, mget(est.names))
    res <- sapply(all.est, comp.crit, signs = signs)
    rownames(res) <- c("sign.recovery", "supp.recovery", "tpp", "fdp", "modelsize", "Nfalse")
    res
}

names(dimnames(results)) <- c("criterion", "estimator", "rep")
dimnames(results)$rep <- 1:R

assign(paste0(type, ".sparsity", s), results)
save(list=paste0(type,".sparsity", s), 
     file=paste0(type,".sparsity", s, ".Rda"))








